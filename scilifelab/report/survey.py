import datetime
import time
import re
import hashlib
import smtplib
from email.mime.text import MIMEText

from scilifelab.db.statusdb import ProjectSummaryConnection
from scilifelab.report.rst import _render, email_templates
from genologics.lims import *
from genologics.config import BASEURI, USERNAME, PASSWORD

def generate_email(email, salt, project):
    """Generate the email text for the survey based on the template
    """
    tpl = email_templates.get('survey')
    hash = hashlib.md5("{}{}".format(salt,project)).hexdigest()
    enddate = (datetime.datetime.now() + datetime.timedelta(days=90)).strftime("%B %d")
    return _render(tpl,**{'hash': hash, 'pname': project})

def send_survey(report, project, email, sender="genomics_support@scilifelab.se", smtphost=None, smtpport=None, dryrun=False):
    """Send the survey to the supplied email address
    """    
    text = generate_email(email, report._meta.salt, project)
    try:
        msg = MIMEText(text, 'html')
        msg['To'] = ",".join(email)
        msg['Subject'] = "Please give feedback about the NGI Sweden service for your project {}".format(project)
        msg['From'] = sender
        if not dryrun:
            s = smtplib.SMTP(host=smtphost, port=smtpport)
            s.sendmail(msg['From'], email, msg.as_string())
            s.quit()  
    except Exception, e:
        report.log.error(e)
        return False
    return True
        
def survey_sent(project):
    """Return True if the field for user survey sent is set, fasle otherwise
    """
    try:
        if type(project.udf.__getitem__('Survey sent')) == datetime.date:
            return True
    except:
        pass
    return False

def project_closed(project, format="%Y-%m-%d"):
    """Return a datetime object representing the close date of the project or None if 
    no information about closing date was found
    """
    date = None
    try:
        date = project.close_date
    except:
        pass
    return date
 
def lims_project(report, pid):
    """Connect to LIMS and get the a project instance
    """
    report.log.debug("Connecting to LIMS instance at {} using username {}".format(BASEURI,USERNAME))
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    if not lims:
        report.log.error("Could not get a connection to LIMS {} API for user {}".format(BASEURI,USERNAME))
        return None
    
    report.log.debug("Getting object for project {}".format(pid))
    project = Project(lims,id=pid)
    return project
    
def project_email(report, project):
    """Obtain the project-associated email addresses from a LIMS project instance
    """
    return [project.researcher.email]
    
def (report, project, **kw):
     
    # Get a connection to the database
    pcon = ProjectSummaryConnection(**kw)
    if not pcon:
        report.log.error("Could not get connection to database".format(project))
        return False
    
    # Get the document for the project
    pdoc = pcon.get_entry(project)
    if not pdoc:
        report.log.error("No such project: {} in database".format(project))
        return False

    # get a project instance from lims
    lproj = lims_project(report, pdoc.get("project_id"))
    if not lproj:
        report.log.error("Could not initiate LIMS object for project {}".format(project))
        return False
        
    # check if project is closed
    closed = project_closed(lproj)
    if closed is None:
        report.log.warn("Project {} is not closed".format(project))
        return False
    report.log.debug("Project {} closed on {}".format(project,datetime.datetime.strptime(closed,report._meta.date_format)))

    # check if a user survey has already been sent
    if survey_sent(lproj):
        report.log.info("Survey already sent for project {}".format(project))
        return False
    report.log.debug("No previous survey sent for {}".format(project))
    
    # get email addresses for persons connected to the project
    emails = project_email(report,lproj)
    if len(emails) == 0:
        report.log.warn("No email addresses found associated with project {}".format(project))
        return False
    
    # verify the format of the email address
    recipients = []
    for email in emails:
        if email is None or not re.match(r'^[_a-z0-9-]+(\.[_a-z0-9-]+)*@[a-z0-9-]+(\.[a-z0-9-]+)*(\.[a-z]{2,3})$',email):
            report.log.warn("Illegal email format: {}".format(email))
            continue
        recipients.append(email)
    
    # send the survey email to each recipient 
    sent = send_survey(report, 
                       project, 
                       recipients, 
                       sender = kw.get("sender"),
                       smtphost=kw.get("smtphost"),
                       smtpport=kw.get("smtpport"),
                       dryrun=report.pargs.dry_run)
    
    # update the project udf to indicate that we have sent out the survey
    if sent:
        report.log.info("Survey sent to recipients {} successfully".format(",".join(recipients)))
        lproj.udf.__setitem__('Survey sent',datetime.datetime.now().date())
        if not report.pargs.dry_run:
            lproj.put()
    elif not sent:
        report.log.warn("Sending survey to recipients {} failed".format(",".join(recipients)))
        
    return sent


        
    
        
    
    