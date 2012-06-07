
import bcbio.log.version
import bcbio.pipeline.config_loader
import sys
import datetime

config = sys.argv[1]
 
versions = bcbio.log.version.get_versions(bcbio.pipeline.config_loader.load_config(config))
versions['bcbb'] = bcbio.log.version._get_git_commit()

print datetime.datetime.now().isoformat()
for sw, ver in versions.items():
    print "%s\t%s" % (sw,ver)
 