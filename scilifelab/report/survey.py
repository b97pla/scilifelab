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

def generate_email(email):
    """Generate the email text for the survey based on the template
    """
    tpl = email_templates.get('survey')
    hash = hashlib.md5("{}{}".format(time.time(),email)).hexdigest()
    enddate = (datetime.datetime.now() + datetime.timedelta(days=90)).strftime("%B %d")
    return _render(tpl,**{'hash': hash, 'enddate': enddate})

def send_survey(report, email, subject="SciLifeLab Survey", sender="genomics_support@scilifelab.se", smtphost="localhost", dryrun=False):
    """Send the survey to the supplied email address
    """
    text = generate_email(email)
    try:
        msg = MIMEText(text)
        msg['To'] = email
        msg['Subject'] = subject
        msg['From'] = sender
        if not dryrun:
            s = smtplib.SMTP(smtphost)
            s.sendmail(msg['From'], recipients, msg.as_string())
            s.quit()  
    except Exception, e:
        report.log.error(e)
        return False
    return True
        
def survey_sent(pdoc):
    """Return True if the field for user survey sent is set to True, fasle otherwise
    """
    return (pdoc.get("user_survey_sent","").lower() == "true")

def project_closed(pdoc, format="%Y-%m-%d"):
    """Return a datetime object representing the close date of the project or None if 
    no information about closing date was found
    """
    date = None
    try:
        date = datetime.datetime.strptime(pdoc.get("close_date"),format)
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
    

    
def initiate_survey(report, project=None, url=None, username=None, password=None):
     
    # Get a connection to the database
    pcon = ProjectSummaryConnection(url=url,username=username,password=password)
    if not pcon:
        report.log.error("Could not get connection to database".format(project))
        return False
    
    # Get the document for the project
    pdoc = pcon.get_entry(project)
    if not pdoc:
        report.log.error("No such project: {} in database".format(project))
        return False
    
    # check if project is closed
    closed = project_closed(pdoc)
    if closed is None:
        report.log.warn("Project {} is not closed".format(project))
        return False
    report.log.debug("Project {} closed on {}".format(project,closed.strftime(report._meta.date_format)))

    # check if a user survey has already been sent
    if survey_sent(pdoc):
        report.log.info("Survey already sent for project {}")
        return False
    report.log.debug("No previous survey sent for {}".format(project))
    
    # get a project instance from lims
    lproj = lims_project(report, pdoc.get("project_id"))
    if not lproj:
        report.log.error("Could not initiate LIMS object for project {}".format(project))
        return False
    
    # get email addresses for persons connected to the project
    emails = project_email(report,lproj)
    if len(emails) == 0:
        report.log.warn("No email addresses found associated with project {}".format(project))
        return False
    
    # send the survey email to each recipient
    succeeded = False
    for email in emails:
        email = 'pontus.larsson@scilifelab.se'
        if email is None or not re.match(r'^[_a-z0-9-]+(\.[_a-z0-9-]+)*@[a-z0-9-]+(\.[a-z0-9-]+)*(\.[a-z]{2,3})$',email):
            report.log.warn("Illegal email format: {}".format(email))
            continue
        
        sent = send_survey(report, email, dryrun=report.pargs.dry_run)
        if not sent:
            report.log.warn("Sending survey to recipient {} failed".format(email))
        succeeded = succeeded or sent
        
    return succeeded

        
    
        
    
    