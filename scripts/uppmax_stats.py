import sys
import datetime
import pexpect

vecAlloc = []
vecUsed = []
vecMem = []

def usage():
    print "python computeStats.py"

def getTimeFromString(string):
    time = string.split(':')
    seconds = 0
    if time[2]!='+':
        seconds = int(time[2])
    seconds += 60*int(time[1])
    if time[0].find('-') != -1:
        #This has DAY information too...
        d = time[0].split('-')
        seconds += 60*60*int(d[1])
        seconds += 24*60*60*int(d[0])
    else:
        seconds += 60*60*int(time[0])
    #dTime = datetime.time(int(time[0]),int(time[1]),int(time[2]))
    return seconds

def main():
    totalAlloc = 0
    totalUsed = 0
    totalMem = 0
    f = open('time','r')
    #Skip header
    readBogus = False
    f.readline()
    f.readline()
    while True:
        if not readBogus:
            l1 = f.readline()
            readBogus = False
        l2 = f.readline()
        if not l1 or not l2:
            break
        c1 = l1.split()
        c2 = l2.split()
        try:
            totalRunTime = getTimeFromString(c2[20])
        except IndexError:
            readBogus = True
            l1 = l2
            continue

        if c1[5] != 'COMPLETED' and c1[5]!= 'TIMEOUT':
            continue

        minProcTime = getTimeFromString(c2[14])
        avgProcTime = getTimeFromString(c2[17])
        allocatedProcs = int(c1[3])
        memory = 0
        try:
            if c2[5].find('G')!=-1:
                memory = float(c2[5].rstrip('G'))
            elif c2[5].find('M')!=-1:
                memory = float(c2[5].rstrip('M'))/1024
            else:
                memory = float(c2[5].rstrip('K'))/(1024*1024)
        except ValueError:
            #This has some weird measure unit. Ignore it.
            print c2[5]
            continue
        if memory > 4.0:
            #Someone needed a fat node, at to themselves
            totalMem += avgProcTime*allocatedProcs
            totalAlloc += avgProcTime*allocatedProcs
            #Don't count these times!
            continue
        if allocatedProcs==1:
            totalAlloc += totalRunTime
            totalUsed += totalRunTime
        else:
            totalAlloc += avgProcTime*allocatedProcs
            if avgProcTime > (totalRunTime)/(allocatedProcs/2):
                #This wasnt's actually parallel
                totalUsed += avgProcTime;
            else:
                totalUsed += avgProcTime*allocatedProcs

    vecAlloc.append(float(totalAlloc)/60)
    vecUsed.append(float(totalUsed)/60)
    vecMem.append(float(totalMem)/60)
    #print 'Total allocated: '+str(totalAlloc/60)+' min'
    #print 'Total used: '+str(totalUsed/60)+' min'

if __name__ == "__main__":

    start = datetime.datetime(2011,1,1,0,0,0)
    step = datetime.timedelta(10)
    end = datetime.datetime(2012,4,1)
    while start < end:
        next = start + step
        sString = start.isoformat('#').split('#')[0]
        eString = next.isoformat('#').split('#')[0]
        runSACCT = 'sacct -S '+sString+' -E '+eString+' -l --alluser -A a2010001,b2010042,b2010045,p2010034,b2010052,b2010065,b2011001,b2011006,b2011011,b2011029,b2011092,b2011163,b2011168,b2011223,b2010062,b2010029,a2010002,a2010003,a2012043 > time '
        rSacct = pexpect.spawn('/bin/bash',['-c', runSACCT])
        rSacct.wait()
        start = next
        main()

    o = open('usageMetric.txt','w')
    for i in range(len(vecAlloc)):
        o.write(str(vecUsed[i])+'\t'+str(vecAlloc[i])+'\t'+str(vecMem[i])+'\n')
    o.close()
