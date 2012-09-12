#!/usr/bin/env python

def main():
	import drmaa
	import os
	import time

	"""Submit a job and wait for it to finish.
	Note, need file called sleeper.sh in home directory.
	"""
	s = drmaa.Session()
	print "Before init:"
	print 'Supported contact strings: ' + s.contact
	print 'Supported DRM systems: ' + str(s.drmsInfo)
	print 'Supported DRMAA implementations: ' + str(s.drmaaImplementation)

	s.initialize()
	print "After init:"
	print 'A DRMAA object was created'
	print 'Supported contact strings: ' + s.contact
	print 'Supported DRM systems: ' + str(s.drmsInfo)
	print 'Supported DRMAA implementations: ' + str(s.drmaaImplementation)
	print 'Version ' + str(s.version)

	jt = s.createJobTemplate()
	jt.jobName = "GalaxyJob"
	jt.remoteCommand = 'echo'
	jt.args = ['foo']
	jt.nativeSpecification = "-A a2010002 -p devel -t 00:00:05"
	jt.workingDirectory = drmaa.JobTemplate.HOME_DIRECTORY
	jt.outputPath = ":"+drmaa.JobTemplate.HOME_DIRECTORY+'/job_stdout.out'
	jt.joinFiles=True # Joins stdout & stderr together
	jt.email=["roman@scilifelab.se"]
	jt.blockEmail=1

	jobid = s.runJob(jt)
	print 'Your job has been submitted with id ' + jobid

	decodestatus = {
		drmaa.JobState.UNDETERMINED: 'process status cannot be determined',
		drmaa.JobState.QUEUED_ACTIVE: 'job is queued and active',
		drmaa.JobState.SYSTEM_ON_HOLD: 'job is queued and in system hold',
		drmaa.JobState.USER_ON_HOLD: 'job is queued and in user hold',
		drmaa.JobState.USER_SYSTEM_ON_HOLD: 'job is queued and in user and system hold',
		drmaa.JobState.RUNNING: 'job is running',
		drmaa.JobState.SYSTEM_SUSPENDED: 'job is system suspended',
		drmaa.JobState.USER_SUSPENDED: 'job is user suspended',
		drmaa.JobState.DONE: 'job finished normally',
		drmaa.JobState.FAILED: 'job finished, but failed',
	}

	status = drmaa.JobState.UNDETERMINED
	while status not in [drmaa.JobState.FAILED, drmaa.JobState.DONE]:
		status = s.jobStatus(jobid)
		print decodestatus[status]
		time.sleep(0.5)

	retval = s.wait(jobid, drmaa.Session.TIMEOUT_NO_WAIT)
	print 'Job: ' + str(retval.jobId) + ' finished with status ' + str(retval.hasExited)
	print """\n
id:                        %(jobId)s
exited:                    %(hasExited)s
signaled:                  %(hasSignal)s
with signal (if signaled): %(terminatedSignal)s
dumped core:               %(hasCoreDump)s
aborted:                   %(wasAborted)s
resource usage:

%(resourceUsage)s
""" % retval._asdict()

	print 'Cleaning up'
	s.deleteJobTemplate(jt)
	s.exit()

	exit(0)
	
if __name__ == "__main__":
	main()
