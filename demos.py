#!/usr/bin/env python
# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#
# author Sergey Lyskov
#

import os, sys, shutil, threading, subprocess, signal, time, re, random, datetime, os.path
import json
from os import path
from optparse import OptionParser, IndentedHelpFormatter


def main(argv):
    ''' Script to run Demos script for Rosetta project '''

    parser = OptionParser(usage="usage: %prog [OPTIONS] [TESTS]", formatter=ParagraphHelpFormatter())
    parser.set_description(main.__doc__)

    parser.add_option("-d", "--database",
      default="", # processed below
      help="Path to Rosetta database. (default: $ROSETTA_DB, ../main/database)",
    )

    parser.add_option("-s", "--source",
      #default=path.join( path.expanduser("~"), "mini"),
      default= path.join( path.dirname( path.dirname(path.abspath(sys.argv[0])) ), 'main/source'),
      help="Directory where Rosetta source repository is (default: ../main/source/)",
    )

    parser.add_option("--tools",
      #default=path.join( path.expanduser("~"), "mini"),
      default= path.join( path.dirname( path.dirname(path.abspath(sys.argv[0])) ), 'tools'),
      help="Directory where Rosetta tools repository is (default: ../tools/)",
    )

    parser.add_option("-j", "--jobs",
      default=1,
      type="int",
      help="number of processors to use on local machine (default: 1)",
    )

    parser.add_option("-t", "--timeout",
      default=0,
      type="int",
      help="Maximum runtime for each test, in minutes (default: no limit)",
      metavar="MINUTES",
    )

    parser.add_option("-c", "--compiler",
      default="gcc",
      help="In selecting binaries, which compiler was used? (default: gcc)",
    )

    parser.add_option("--mode",
      default="release",
      help="In selecting binaries, which mode was used? (default: release)",
    )

    parser.add_option("--extras",
      default="default",
      dest="extras",
      help="in selecting binaries, which options were specified? (default: default)",
    )

    # parser.add_option("--daemon", action="store_true", dest="daemon", default=False,
    #   help="generate daemon friendly output (off by default)"
    # )

    parser.add_option("--additional_flags",
      default="",
      help="Add additional flags to integration tests. (default: None)",
    )

    (options, remaining_args) = parser.parse_args(args=argv)
    global Options;  Options = options
    Options.num_procs = Options.jobs

    # Strip off whitespace and blank remaining arguments
    args = [a.strip() for a in remaining_args if a.strip()]

    options.source = path.abspath( options.source )
    print 'Using Rosetta source dir at:', options.source

    options.tools = path.abspath( options.tools )
    print 'Using Rosetta tools dir at:', options.tools

    if options.database == parser.get_default_values().database:
        if os.environ.get('ROSETTA3_DB') is not None and \
                path.isdir(os.environ.get('ROSETTA3_DB')):
            options.database = os.environ.get('ROSETTA3_DB')
        else:
            options.database = path.join( path.dirname( path.dirname(path.abspath(sys.argv[0])) ), 'main/database')
            print options.database

        if not path.isdir( options.database ):
            options.database = path.join( path.expanduser("~"), "rosetta_database")

        if not path.isdir( options.database ):
            print "Can't find database at %s; please set $ROSETTA3_DB or use -d" % options.database
            return 1
    # Normalize path before we change directories!
    options.database = path.abspath(options.database)

    # Make sure the current directory is the script directory:
    # Using argv[] here causes problems when people try to run the script as "python integration.py ..."
    #os.chdir( path.dirname(sys.argv[0]) ) # argv[0] is the script name
    if not path.isdir("public"):
        print "You must run this script from rosetta/demos/"
        return 2


    outdir = "public.run"
    if path.isdir(outdir): print 'Removing old run-dir: %s...' % outdir;  shutil.rmtree(outdir)  # remove old dir if any
    os.mkdir(outdir)


    if len(args) > 0: tests = args
    else: tests = [ d for d in os.listdir("public") if not d.startswith(".") and path.isdir(path.join("public", d)) ]

    def signal_handler(signal_, f):
        print 'Ctrl-C pressed... killing child jobs...'
        for nt in Jobs:
            os.killpg(os.getpgid(nt.pid), signal.SIGKILL)

    signal.signal(signal.SIGINT, signal_handler)


    runtimes={}
    queue = Queue()
    queue.TotalNumberOfTasks = len(tests)

    # Write substitution parameters to result directory
    with open(path.join( outdir, "test_parameters.json"), "w") as parameters_file:
        json.dump(generateIntegrationTestGlobalSubstitutionParameters(), parameters_file, sort_keys=True, indent=2)

    for test in tests:
        #commandfile = path.join("public",test,"/command")
        commandfilepath = path.join("public",test) + '/command'

        if(os.path.isfile(commandfilepath)):
            print '### This exists: '+commandfilepath  ## LGN 20131108
            queue.put(test)
            #shutil.copytree( path.join("tests", test), path.join(outdir, test) )
            print '~~~', path.join("public", test), path.join(outdir, test)
        
            copytree( path.join("public", test), path.join(outdir, test) )  #  accept=lambda src, dst: path.basename(src) != '.svn' )

    while not queue.empty():
        test = queue.get()
        if test is None: break

        cmd_line_sh, workdir = generateIntegrationTestCommandline(test, outdir);

        def run(times):
            #execute('Running Test %s' % test, 'bash ' + cmd_line_sh)
            extra = 'ulimit -t%s && ' % Options.timeout  if Options.timeout else ''
            res = execute('Running Test %s' % test, '%sbash %s' % (extra, cmd_line_sh), return_=True)
            if res:
                error_string = "*** Test %s did not run!  Check your --mode flag and paths. [%s]\n" % (test, datetime.datetime.now())
                file(path.join(nt.workdir, ".test_did_not_run.log"), 'w').write(error_string)
                print error_string,
                times[test] = float('nan')

        def normal_finish(nt, times):
            queue.task_done()
            percent = (100* (queue.TotalNumberOfTasks-queue.qsize())) / queue.TotalNumberOfTasks
            elapse_time = time.time() - nt.start_time
            print "Finished %-40s in %3i seconds\t [~%4s test (%s%%) started, %4s in queue, %4d running]" % (nt.test, elapse_time, queue.TotalNumberOfTasks-queue.qsize(), percent, queue.qsize(), queue.unfinished_tasks-queue.qsize() )
            if nt.test not in times:
                times[nt.test] = elapse_time

        def error_finish(nt, times):
            error_string = "*** Test %s did not run!  Check your --mode flag and paths. [%s]\n" % (test, datetime.datetime.now())
            file(path.join(nt.workdir, ".test_did_not_run.log"), 'w').write(error_string)
            print error_string,
            times[nt.test] = float('nan')
            normal_finish(nt, times)

        def timeout_finish(nt, times):
            error_string = "*** Test %s exceeded the timeout=%s  and will be killed! [%s]\n" % (test, Options.timeout, datetime.datetime.now())
            file(path.join(nt.workdir, ".test_got_timeout_kill.log"), 'w').write(error_string)
            print error_string,
            times[nt.test] = float('inf')
            normal_finish(nt, times)

        if Options.jobs > 1:
            pid, nt = mFork(times=runtimes, test=test, workdir=workdir, queue=queue, timeout=Options.timeout, normal_finish=normal_finish, error_finish=error_finish, timeout_finish=timeout_finish)
            if not pid:  # we are child process
                signal.signal(signal.SIGINT, signal.SIG_DFL)
                run(runtimes)
                sys.exit(0)
        else:
            nt = NT(times=runtimes, test=test, workdir=workdir, queue=queue, start_time=time.time(), timeout=Options.timeout, normal_finish=normal_finish, error_finish=error_finish, timeout_finish=timeout_finish)
            run(runtimes)
            if nt.timeout and (time.time() - nt.start_time > nt.timeout): nt.timeout_finish(nt)
            else: normal_finish(nt,runtimes)

    mWait(all_=True)  # waiting for all jobs to finish before movinf in to next phase


# -------------------------------------
class NT:  # named tuple
    def __init__(self, **entries): self.__dict__.update(entries)
    def __repr__(self):
        r = '|'
        for i in dir(self):
            if not i.startswith('__'): r += '%s --> %s, ' % (i, getattr(self, i))
        return r[:-2]+'|'

Jobs = []  # Global list of NameTuples  (pid, tag, start_time, out_dir,...)

def write_runtimes(runtimes, dir):
    try:
      time_file = open(dir+'/runtimes.yaml', 'w')
      json.dump(runtimes, time_file, sort_keys=True, indent=2)
      time_file.close()
    except:
      pass # if no JSON, just forget this step



def execute(message, command_line, return_=False, untilSuccesses=False, print_output=True, verbose=True):
    if verbose:
        print message
        print command_line

    while True:
        #(res, output) = commands.getstatusoutput(commandline)

        po = subprocess.Popen(command_line+ ' 1>&2', bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #po = subprocess.Popen(command_line+ ' 1>&2', bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        f = po.stderr
        output = ''
        for line in f:
            #po.poll()
            if print_output: print line,
            output += line
            sys.stdout.flush()
        f.close()
        while po.returncode is None: po.wait()
        res = po.returncode
        #print '_____________________', res

        if res and untilSuccesses: pass  # Thats right - redability COUNT!
        else: break

        print "Error while executing %s: %s\n" % (message, output)
        print "Sleeping 60s... then I will retry..."
        time.sleep(60)

    if res:
        if print_output: print "\nEncounter error while executing: " + command_line
        if not return_: sys.exit(1)

    if return_ == 'output': return output
    else: return res

def print_(msg, color=None, background=None, bright=False, blink=False, action='print', endline=True):
    ''' print string with color and background. Avoid printing and return results str instead if action is 'return'. Also check for 'Options.no_color'
    '''
    colors = dict(black=0, red=1, green=2, yellow=3, blue=4, magenta=5, cyan=6, white=7)  # standard ASCII colors

    if 'Options' in globals()  and  hasattr(Options, 'color')  and  not Options.color: s = str(msg)
    else:
        s  = ['3%s' % colors[color] ] if color else []
        s += ['4%s' % colors[background] ] if background else []
        s += ['1'] if bright else []
        s += ['5'] if blink else []

        s = '\033[%sm%s\033[0m%s' % (';'.join(s), msg, '\n' if endline else '')

    if action == 'print': sys.stdout.write(s)
    else: return s


def mFork(times, tag=None, overhead=0, **args):
    ''' Check if number of child process is below Options.jobs. And if it is - fork the new pocees and return its pid.
    '''
    #print_('Groups:%s' % os.getgroups(), color='cyan')
    while len(Jobs) >= Options.jobs + overhead:
        for j in Jobs[:] :
            r = os.waitpid(j.pid, os.WNOHANG)
            if r == (j.pid, 0):  # process have ended without error
                normal_finish = getattr(j, 'normal_finish', lambda x: None)
                normal_finish(j, times)
                Jobs.remove(j)

            elif r[0] == j.pid :
                error_finish = getattr(j, 'error_finish', lambda x: None)
                error_finish(j, times)
                Jobs.remove(j)

            else:
                #pass
                if j.timeout:
                    if time.time() - j.start_time > j.timeout :
                        #print '~~~~~~~~~~~ pids:', j.pid, os.getpid(), os.getppid()
                        #print '~~~~~~~~~~~ groups:', os.getpgid(j.pid), os.getpgrp()
                        os.kill(j.pid, signal.SIGKILL)  #
                        #os.killpg(os.getpgid(j.pid), signal.SIGKILL)
                        timeout_finish = getattr(j, 'timeout_finish', lambda x: None)
                        timeout_finish(j, times)
                        Jobs.remove(j)
                        break

        if len(Jobs) >= Options.jobs + overhead: time.sleep(.2)

    sys.stdout.flush();  sys.stderr.flush();
    pid = os.fork()
    if pid: pass # We are parent!
    Jobs.append( NT(times=times, pid=pid, tag=tag, start_time=time.time(), **args) )
    return pid, Jobs[-1]


def mWait(tag=None, all_=False, timeout=0):
    ''' Wait for process tagged with 'tag' for completion
    '''
    while True :
        #print 'Waiting for %s: ' % tag, Jobs

        for j in [ x for x in Jobs if x.tag==tag or all_==True]:
            #print 'Waiting2: ', Jobs
            #try:
            r = os.waitpid(j.pid, os.WNOHANG)
            times = j.times
            if r == (j.pid, 0):  # process have ended without error
                normal_finish = getattr(j, 'normal_finish', lambda x: None)
                normal_finish(j, times)
                Jobs.remove(j)

            elif r[0] == j.pid :  # process ended but with error, special case we will have to wait for all process to terminate and call system exit.
                error_finish = getattr(j, 'error_finish', lambda x: None)
                error_finish(j, times)
                Jobs.remove(j)

            elif j.timeout:
                if time.time() - j.start_time > j.timeout :
                    os.kill(j.pid, signal.SIGKILL)  #os.killpg(os.getpgid(j.pid), signal.SIGKILL)
                    timeout_finish = getattr(j, 'timeout_finish', lambda x: None)
                    timeout_finish(j, times)
                    Jobs.remove(j)

        time.sleep(.2)
        if not Jobs: return


def generateIntegrationTestGlobalSubstitutionParameters(host=None):
    # Variables that may be referrenced in the cmd string:
    python = sys.executable
    source = Options.source
    database = Options.database
    tools = Options.tools

    bin = path.join(source, "bin")
    pyapps = path.join(source, "src", "python", "apps")

    if sys.platform.startswith("linux"):
        platform = "linux" # can be linux1, linux2, etc
    elif sys.platform == "darwin":
        platform = "macos"
    elif sys.platform == "cygwin":
        platform = "cygwin"
    else:
        platform = "_unknown_"

    compiler = Options.compiler
    mode = Options.mode
    extras = Options.extras
    binext = Options.extras+"."+platform+compiler+mode
    additional_flags = Options.additional_flags
    # dbms_host = Options.dbms_host
    # dbms_user = Options.dbms_user
    # dbms_port = Options.dbms_port

    return dict(locals())

def generateIntegrationTestSubstitutionParameters(test, outdir, host=None):
    """ Generate substitution parameters for integration command generation."""

    params = generateIntegrationTestGlobalSubstitutionParameters(host)
    # params["dbms_database_name"] = Options.dbms_database_name % { 'test': test }
    # params["dbms_pq_schema"] = Options.dbms_pq_schema % { 'test': test }
    params["workdir"] = path.abspath( path.join(outdir, test) )

    return params

def generateIntegrationTestCommandline(test, outdir, host=None):
    ''' Generate and write command.sh and return command line that will run given integration test
    '''
    # Read the command from the file "command"
    params = generateIntegrationTestSubstitutionParameters(test, outdir, host)
    workdir = params["workdir"]

    cmd=''
    # A horrible hack b/c SSH doesn't honor login scripts like .bash_profile
    # when executing specific remote commands.
    # This causes problems with e.g. the custom Python install on the Whips.
    # So we replace the default remote PATH with the current local one.
    if host is not None:
      cmd = 'PATH="%s"\n%s' % (os.environ["PATH"], cmd)
    cmd += '\n'
    cmd += file(path.join(workdir, "command")).read().strip()
    cmd = cmd % params # variable substitution using Python printf style

    cmd_line_sh = path.join(workdir, "command.sh")
    f = file(cmd_line_sh, 'w');  f.write(cmd);  f.close() # writing back so test can be easily re-run by user lately...
    #if "'" in cmd: raise ValueError("Can't use single quotes in command strings!")
    #print cmd; print

    return cmd_line_sh, workdir

class Worker:
    def __init__(self, queue, outdir, opts, times, host=None, timeout_minutes=0):
        self.queue = queue
        self.outdir = outdir
        self.opts = opts
        self.host = host
        self.timeout = timeout_minutes * 60
        self.times = times

    def work(self):
        running=0
        try:
            while True:
                test = self.queue.get_nowait()
                try: # Actually catch exception and ignore it.  Python 2.4 can't use "except" and "finally" together.
                    start = time.time() # initial guess at start time, in case of exception
                    try: # Make sure job is marked done even if we throw an exception
                        cmd_line_sh, workdir = generateIntegrationTestCommandline(test, self.outdir, host=self.host)

                        if self.host is None:
                            print "Running  %-40s on localhost ..." % test
                            proc = subprocess.Popen(["bash",  cmd_line_sh], preexec_fn=os.setpgrp)
                        # Can't use cwd=workdir b/c it modifies *local* dir, not remote dir.
                        else:
                            print "Running  %-40s on %20s ..." % (test, self.host)
                            bash_cmd='bash '+cmd_line_sh
                            proc = subprocess.Popen(["ssh", self.host, bash_cmd], preexec_fn=os.setpgrp)#, cwd=workdir)
                            #os._exit(os.EX_IOERR)
                        start = time.time() # refined start time
                        if self.timeout == 0:
                            retcode = proc.wait() # does this block all threads?
                        else:
                            while time.time() - start <= self.timeout:
                                retcode = proc.poll()
                                if retcode is not None: break
                                time.sleep(1)
                            if retcode is None:
                                print "*** Test %s exceeded the timeout and will be killed! [%s]\n" % (test, datetime.datetime.now())
                                self.times[test] = float('inf')
                                #os.kill(proc.pid, signal.SIGTERM)
                                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                        if retcode != 0 and retcode is not None:
                            self.times[test] = float('nan')
                            if self.host is None:
                              error_string = "*** Test %s did not run on host %s!  Check your --mode flag and paths. [%s]\n" % (test, 'local_host', datetime.datetime.now())
                            else:
                              error_string = "*** Test %s did not run on host %s!  Check your --mode flag and paths. [%s]\n" % (test, self.host, datetime.datetime.now())
                            print error_string,

                            # Writing error_string to a file, so integration test should fail for sure
                            file(path.join(workdir, ".test_did_not_run.log"), 'w').write(error_string)

                    finally: # inner try
                        percent = (100* (self.queue.TotalNumberOfTasks-self.queue.qsize())) / self.queue.TotalNumberOfTasks
                        elapse_time = time.time() - start
                        print "Finished %-40s in %3i seconds\t [~%4s test (%s%%) started, %4s in queue, %4d running]" % (test, elapse_time, self.queue.TotalNumberOfTasks-self.queue.qsize(), percent, self.queue.qsize(), self.queue.unfinished_tasks-self.queue.qsize() )
                        if test not in self.times: self.times[test] = elapse_time
                        self.queue.task_done()

                except Exception, e: # middle try
                    print e
        except Empty: # outer try
            pass # we're done, just return


class ParagraphHelpFormatter(IndentedHelpFormatter):
    '''
    A help formatter that respects paragraph breaks (blank lines) in usage strings.
    '''
    def _format_text(self, text):
        paragraphs = re.split('\n([ \t]*\n)+', text)
        paragraphs = [ IndentedHelpFormatter._format_text(self, p.strip()) for p in paragraphs ]
        return '\n'.join(paragraphs) # each already ends in a newline

def copytree(src, dst, symlinks=False, accept=lambda srcname, dstname: True):
    """Recursively copy a directory tree using copy2(), with filtering.
    Copied from shutil so I could filter out .svn entries.
    """
    names = os.listdir(src)
    os.makedirs(dst)
    errors = []
    for name in names:
        srcname = os.path.join(src, name)
        dstname = os.path.join(dst, name)
        if not accept(srcname, dstname): continue
        try:
            if symlinks and os.path.islink(srcname):
                linkto = os.readlink(srcname)
                os.symlink(linkto, dstname)
            elif os.path.isdir(srcname):
                copytree(srcname, dstname, symlinks, accept)
            else:
                shutil.copy2(srcname, dstname)
            # XXX What about devices, sockets etc.?
        except (IOError, os.error), why:
            errors.append((srcname, dstname, str(why)))
        # catch the Error from the recursive copytree so that we can
        # continue with other files
        except shutil.Error, err:
            errors.extend(err.args[0])
    try:
        shutil.copystat(src, dst)
    #except WindowsError:
        # can't copy file access times on Windows
        pass
    except OSError, why:
        errors.extend((src, dst, str(why)))
    if errors:
        raise shutil.Error, errors

################################################################################
# Python 2.4 lacks support for join() / task_done() in the Queue class,
# so I pasted the 2.5 implementation here.
# With 2.5+, you can just do "from Queue import *" instead.

from time import time as _time
from collections import deque

class Empty(Exception):
    "Exception raised by Queue.get(block=0)/get_nowait()."
    pass

class Full(Exception):
    "Exception raised by Queue.put(block=0)/put_nowait()."
    pass

class Queue:
    """Create a queue object with a given maximum size.

    If maxsize is <= 0, the queue size is infinite.
    """
    def __init__(self, maxsize=0):
        try:
            import threading
        except ImportError:
            import dummy_threading as threading
        self._init(maxsize)
        # mutex must be held whenever the queue is mutating.  All methods
        # that acquire mutex must release it before returning.  mutex
        # is shared between the three conditions, so acquiring and
        # releasing the conditions also acquires and releases mutex.
        self.mutex = threading.Lock()
        # Notify not_empty whenever an item is added to the queue; a
        # thread waiting to get is notified then.
        self.not_empty = threading.Condition(self.mutex)
        # Notify not_full whenever an item is removed from the queue;
        # a thread waiting to put is notified then.
        self.not_full = threading.Condition(self.mutex)
        # Notify all_tasks_done whenever the number of unfinished tasks
        # drops to zero; thread waiting to join() is notified to resume
        self.all_tasks_done = threading.Condition(self.mutex)
        self.unfinished_tasks = 0

    def task_done(self):
        """Indicate that a formerly enqueued task is complete.

        Used by Queue consumer threads.  For each get() used to fetch a task,
        a subsequent call to task_done() tells the queue that the processing
        on the task is complete.

        If a join() is currently blocking, it will resume when all items
        have been processed (meaning that a task_done() call was received
        for every item that had been put() into the queue).

        Raises a ValueError if called more times than there were items
        placed in the queue.
        """
        self.all_tasks_done.acquire()
        try:
            unfinished = self.unfinished_tasks - 1
            if unfinished <= 0:
                if unfinished < 0:
                    raise ValueError('task_done() called too many times')
                self.all_tasks_done.notifyAll()
            self.unfinished_tasks = unfinished
        finally:
            self.all_tasks_done.release()

    def join(self):
        """Blocks until all items in the Queue have been gotten and processed.

        The count of unfinished tasks goes up whenever an item is added to the
        queue. The count goes down whenever a consumer thread calls task_done()
        to indicate the item was retrieved and all work on it is complete.

        When the count of unfinished tasks drops to zero, join() unblocks.
        """
        self.all_tasks_done.acquire()
        try:
            while self.unfinished_tasks:
                self.all_tasks_done.wait()
        finally:
            self.all_tasks_done.release()

    def qsize(self):
        """Return the approximate size of the queue (not reliable!)."""
        self.mutex.acquire()
        n = self._qsize()
        self.mutex.release()
        return n

    def empty(self):
        """Return True if the queue is empty, False otherwise (not reliable!)."""
        self.mutex.acquire()
        n = self._empty()
        self.mutex.release()
        return n

    def full(self):
        """Return True if the queue is full, False otherwise (not reliable!)."""
        self.mutex.acquire()
        n = self._full()
        self.mutex.release()
        return n

    def put(self, item, block=True, timeout=None):
        """Put an item into the queue.

        If optional args 'block' is true and 'timeout' is None (the default),
        block if necessary until a free slot is available. If 'timeout' is
        a positive number, it blocks at most 'timeout' seconds and raises
        the Full exception if no free slot was available within that time.
        Otherwise ('block' is false), put an item on the queue if a free slot
        is immediately available, else raise the Full exception ('timeout'
        is ignored in that case).
        """
        self.not_full.acquire()
        try:
            if not block:
                if self._full():
                    raise Full
            elif timeout is None:
                while self._full():
                    self.not_full.wait()
            else:
                if timeout < 0:
                    raise ValueError("'timeout' must be a positive number")
                endtime = _time() + timeout
                while self._full():
                    remaining = endtime - _time()
                    if remaining <= 0.0:
                        raise Full
                    self.not_full.wait(remaining)
            self._put(item)
            self.unfinished_tasks += 1
            self.not_empty.notify()
        finally:
            self.not_full.release()

    def put_nowait(self, item):
        """Put an item into the queue without blocking.

        Only enqueue the item if a free slot is immediately available.
        Otherwise raise the Full exception.
        """
        return self.put(item, False)

    def get(self, block=True, timeout=None):
        """Remove and return an item from the queue.

        If optional args 'block' is true and 'timeout' is None (the default),
        block if necessary until an item is available. If 'timeout' is
        a positive number, it blocks at most 'timeout' seconds and raises
        the Empty exception if no item was available within that time.
        Otherwise ('block' is false), return an item if one is immediately
        available, else raise the Empty exception ('timeout' is ignored
        in that case).
        """
        self.not_empty.acquire()
        try:
            if not block:
                if self._empty():
                    raise Empty
            elif timeout is None:
                while self._empty():
                    self.not_empty.wait()
            else:
                if timeout < 0:
                    raise ValueError("'timeout' must be a positive number")
                endtime = _time() + timeout
                while self._empty():
                    remaining = endtime - _time()
                    if remaining <= 0.0:
                        raise Empty
                    self.not_empty.wait(remaining)
            item = self._get()
            self.not_full.notify()
            return item
        finally:
            self.not_empty.release()

    def get_nowait(self):
        """Remove and return an item from the queue without blocking.

        Only get an item if one is immediately available. Otherwise
        raise the Empty exception.
        """
        return self.get(False)

    # Override these methods to implement other queue organizations
    # (e.g. stack or priority queue).
    # These will only be called with appropriate locks held

    # Initialize the queue representation
    def _init(self, maxsize):
        self.maxsize = maxsize
        self.queue = deque()

    def _qsize(self):
        return len(self.queue)

    # Check whether the queue is empty
    def _empty(self):
        return not self.queue

    # Check whether the queue is full
    def _full(self):
        return self.maxsize > 0 and len(self.queue) == self.maxsize

    # Put a new item in the queue
    def _put(self, item):
        self.queue.append(item)

    # Get an item from the queue
    def _get(self):
        return self.queue.popleft()
################################################################################


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
