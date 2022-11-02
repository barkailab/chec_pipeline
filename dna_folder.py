import configparser, os, sys, glob, time, argparse, getpass
from datetime import datetime
## get user input
parser = argparse.ArgumentParser(description='settings for the file selection')
parser.add_argument('--infolder', default=os.getcwd())
parser.add_argument('--outfolder', default='')
parser.add_argument('--tempfolder', default='')
parser.add_argument('--config', default='/home/labs/barkailab/LAB/scripts/felix_bcl2fastq/std_opts.txt')
parser.add_argument('--ending', default='_R1_001.fastq.gz')
parser.add_argument('--cores', default='4')
parser.add_argument('--queue', default='molgen-q')
parser.add_argument('--memSize', default='8000')
parser.add_argument('--overWrite',type=bool, default=True)

args = parser.parse_args()
##generate run specific config file
runid=getpass.getuser()+datetime.now().strftime("%y%m%d");
config = configparser.ConfigParser() 
if '/' in args.config:
	config.read(args.config)
else:
	config.read('/home/labs/barkailab/LAB/scripts/felix_bcl2fastq/'+args.config)

	
if not config.has_section('fileopts'):
	config.add_section('fileopts')
config.set('fileopts','infolder',args.infolder);
if (len(args.outfolder)) > 0:
	config.set('fileopts','outfolder',args.outfolder);
else:
	config.set('fileopts','outfolder',args.infolder);
config.set('fileopts','ending',args.ending);
if (len(args.tempfolder)) > 0:
	config.set('fileopts','tempfolder',args.tempfolder);
else:
	config.set('fileopts','tempfolder',config['fileopts']['outfolder']+'/temp/');
if not os.path.isdir(config['fileopts']['outfolder']):
	os.mkdir(config['fileopts']['outfolder'])
if not os.path.isdir(config['fileopts']['tempfolder']):
	os.mkdir(config['fileopts']['tempfolder'])

config.set('fileopts','statfile',config['fileopts']['outfolder']+'/'+runid+'_stats.tsv')

configfile=config['fileopts']['outfolder']+'/'+runid+'_opts.txt'
if config.has_section('prepro'):
	os.system('echo "File\tfailed filter\tpassed filter\tperfect al.\tmultiple al.\tdisc. al.\tsingle al.\tno al.\t al. pairs\t duplicates\tpairse used\tlibrary size\tmed. fragment\tTarget">'+config['fileopts']['statfile'])
else:
	os.system('echo "File\ttotal reads\tperfect al.\tmultiple al.\tdisc. al.\tsingle al.\tno al.\t al. pairs\t duplicates\tpairs used\tlibrary size\tmed. fragment">'+config['fileopts']['statfile'])

config_write=open(configfile,'w')
config.write(config_write)
config_write.close()


##look up all files and run pipeline
import sys,os,glob
print(config['fileopts']['infolder'])
filenames = glob.glob(config['fileopts']['infolder']+'/*'+config['fileopts']['ending']) # find all files
cmd_ind='bsub -oo ' +config['fileopts']['outfolder']+ '/dna_folder.log -n '+args.cores +' -q ' +args.queue +' -R "span[hosts=1]" -R "rusage[mem=%s]" "python /home/labs/barkailab/LAB/scripts/felix_bcl2fastq/dna_ind.py %s %s"'
cmdNoBsub='python /home/labs/barkailab/LAB/scripts/felix_bcl2fastq/dna_ind.py %s %s';
for i in range(0,len(filenames)):	
	filebase=os.path.basename(filenames[i])	
	skipFile=False
	if not args.overWrite:
		if (os.path.isfile(config['fileopts']['outfolder']+'/'+filebase[:-len(config['fileopts']['ending'])]+'.out')):
			skipFile=True
	if not 'Undetermined' in filebase:
		if not skipFile:
			print(filebase[:-len(config['fileopts']['ending'])])
			if args.queue=='none':
				os.system(cmdNoBsub%(filebase[:-len(config['fileopts']['ending'])],configfile))
			elif os.path.getsize(filenames[i])>400000000:
				os.system(cmd_ind%('16000',filebase[:-len(config['fileopts']['ending'])],configfile))
			else:
				os.system(cmd_ind%(args.memSize,filebase[:-len(config['fileopts']['ending'])],configfile))
print(len(filenames))