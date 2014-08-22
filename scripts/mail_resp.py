
# -*- coding: utf-8 -*-
from genologics.lims import *
from genologics.config import BASEURI, USERNAME, PASSWORD
from datetime import * 
import smtplib
from email.mime.text import MIMEText


lims = Lims(BASEURI, USERNAME, PASSWORD)

sixMonthsAgo=date.today()-timedelta(weeks=26);
yesterday=date.today()-timedelta(days=1)
pjs=lims.get_projects(open_date=sixMonthsAgo.strftime("%Y-%m-%d"))

summary={}
email={u'Francesco Vezzi':'francesco.vezzi@scilifelab.se',
        u'Chuan Wang':'chuan.wang@scilifelab.se',
        u'Remi-Andre Olsen':'remi-andre.olsen@scilifelab.se',
        u'Britta Lotstedt':'britta.lotstedt@scilifelab.se',
        u'Joel Gruselius':'joel.gruselius@scilifelab.se',
        u'Par Lundin':'par.lundin@scilifelab.se',
        u'Nemanja Rilakovic':'nemanja.rilakovic@scilifelab.se',
        u'Senthilkumar Paneerselvam':'senthilkumar.paneerselvam@scilifelab.se',
        u'Carolina Bonilla':'carolina.bonilla@scilifelab.se',
        u'Helena Samuelsson':'helena.samuelsson@scilifelab.se',
        u'Mario Giovacchini':'mario.giovacchini@scilifelab.se',
        u'Marianna Sjogren':'marianna.sjogren@scilifelab.se',
        u'Jun Wang':'jun.wang@scilifelab.se',
        u'Bahram Amini':'bahram.amini@scilifelab.se',
        u'Christian Natanaelsson':'christian.natanaelsson@scilifelab.se',
        u'Mattias Oskarsson':'mattias.oskarsson@scilifelab.se'
        }
project_types=['Illumina Sequencing (Illumina SBS) 4.0', 'MiSeq Run (MiSeq) 4.0','Cluster Generation (Illumina SBS) 4.0','Denature, Dilute and Load Sample (MiSeq) 4.0', 'Aggregate QC (DNA) 4.0','Aggregate QC (RNA) 4.0', 'Project Summary 1.3']

def clean_names(name):
    return name.replace(u"\u00f6", "o").replace(u"\u00e9", "e").replace(u"\u00e4", "a")

for p in pjs:
    #Assuming this will be run on the early morning, this grabs all processes from the list that have been modified the day before
    pro=lims.get_processes(projectname=p.name, type=project_types, last_modified=yesterday.strftime("%Y-%m-%dT00:00:00Z"));
    completed=[]
    bfr=None
    lbr=None
    if pro:
        for pr in pro:
            date_start=None
            #Special case for the project summary
            if pr.type.name== 'Project Summary 1.3': 
                if 'Queued' in pr.udf and pr.udf['Queued'] == yesterday.strftime("%Y-%m-%d"):
                    completed.append({'project':p.name, 
                        'action':'has been queued', 
                        'date':pr.udf['Queued'], 
                        'techID':pr.udf['Signature Queued'],
                        'tech':pr.technician.first_name+" "+pr.technician.last_name, 
                        'sum':True})
                if 'All samples sequenced' in pr.udf and pr.udf['All samples sequenced'] == yesterday.strftime("%Y-%m-%d"):
                    completed.append({'project':p.name, 
                        'action':'Has all its samples sequenced',
                        'date':pr.udf['All samples sequenced'], 
                        'techID':pr.udf['Signature All samples sequenced'],
                        'tech':pr.technician.first_name+" "+pr.technician.last_name, 
                        'sum':True})
                if ' All raw data delivered' in pr.udf and pr.udf[' All raw data delivered'] == yesterday.strftime("%Y-%m-%d"):
                    completed.append({'project':p.name, 
                        'action':'Has all its samples sequenced',
                        'date':pr.udf[' All raw data delivered'], 
                        'techID':pr.udf['Signature  All raw data delivered'],
                        'tech':pr.technician.first_name+" "+pr.technician.last_name, 
                        'sum':True})

            else:#I don't want to combine this in a single elif because of the line 80, that must be done in the else, but regardless of the if
                if 'Run ID' in pr.udf: #this is true for sequencing processes, and gives the actual starting date
                    date_start=pr.udf['Run ID'].split('_')[0]#format is YYMMDD
                    date_start=date_start[:2]+"-"+date_start[2:4]+"-"+date_start[4:6]
                    if date_start==pr.date_run[4:]:
                        date_start=None
                    else:
                        date_start="20"+date_start#now, the format is YYYY-MM-DD, assuming no prjects come from the 1990's or the next century...
                completed.append({'project':p.name, 'process':pr.type.name, 'limsid':pr.id, 'start':date_start, 'end': pr.date_run, 'tech':pr.technician.first_name+" "+pr.technician.last_name,'sum':False}) 
        if len(completed)>0:#If we actually have stuff to mail
            ps=lims.get_processes(projectname=p.name, type='Project Summary 1.3')
            for oneps in ps:#there should be only one project summary per project anyway.
                if 'Bioinfo responsible' in oneps.udf  :
                    bfr=clean_names(oneps.udf['Bioinfo responsible'])

                    summary[bfr]=completed
                if 'Lab responsible' in oneps.udf:
                    lbr=clean_names(oneps.udf['Lab responsible'])
                    summary[lbr]=completed

for resp in summary:
    plist=set()#no duplicates
    body=''
    for struct in summary[resp]:
        if resp != struct['tech'] and not struct['sum']:    
            plist.add(struct['project'])
            body+="In project {},  {} ({})".format(struct['project'], struct['process'], struct['limsid'])
            if struct['start']:
                body+="started on {}, ".format(struct['start'])
            elif struct['end']:
                body+="ended on {}, ".format(struct['end'])
            else:  
                body+="has been modified yesterday, "
            body+="Done by {}\n".format(struct['tech'])
        elif struct['sum']:
            plist.add(struct['project'])
            body+='Project {} {} on {} by {}\n'.format(struct['project'], struct['action'],struct['date'],struct['techID'])
    #print "TO " +email[resp]
    #print " ".join(plist)
    #print body
    if body!= '':
        msg=MIMEText(body)
        msg['Subject']='[Lims update] {}'.format(" ".join(plist))
        msg['From']='Lims_monitor'
        msg['To'] =email[resp]
        s = smtplib.SMTP('smtp.ki.se')
        print msg
        #s.sendmail('genologics-lims@scilifelab.se', msg['To'], msg.as_string())
        s.quit()
