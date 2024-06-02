import configparser,os,subprocess,re,sys,numpy
# this si the pipeline processing individual samples - for better troubleshooting each commadn is first printed before it is actually run using subprocess.check_output() 
config = configparser.ConfigParser()
config.read(sys.argv[2])
infile = config['fileopts']['infolder']+'/'+sys.argv[1]
# The pipeline first checks if alignment and duplciated finding (the most ressource intensive steps) were already performed previously - IF NOT it goes into it. 
if not (os.path.isfile(config['fileopts']['tempfolder']+'/'+sys.argv[1] + '_wd.bam') or os.path.isfile(config['fileopts']['tempfolder']+'/'+sys.argv[1] + 'F_wd.bam')):
	# the pipeline checks if its nececessary to preprocess/filter the data and if yes which program should be used cutadapt for ChEC-seq adapter dimers (TAT) or for ATAC adapters (ATAC)
	if config.has_section('prepro'): 
		prepro=0;
		if 'lc' in config['prepro']['filter']:
			cmd ="/home/labs/barkailab/LAB/scripts/felix_bcl2fastq/prinseq++ -fastq %s_R1_001.fastq.gz -fastq2 %s_R2_001.fastq.gz -"+config['prepro']['filter']+" -out_good %sF_R1_001.fastq -out_good2 %sF_R2_001.fastq -out_single /dev/null/ -out_single2 /dev/null/ -out_bad /dev/null/ -out_bad2 /dev/null/"
			print(cmd %(infile,infile,infile,infile))
			psout=subprocess.check_output(cmd %(infile,infile,infile,infile),shell = True, stderr=subprocess.STDOUT)
			if (type(psout)==bytes):
				psout=psout.decode('utf-8');
			psout=list(map(int,re.findall('\d+',psout)))
			prepro=prepro+psout[0];
			subprocess.check_output("gzip -f %sF_R*_001.fastq" %(infile),shell = True, stderr=subprocess.STDOUT)
			sys.argv[1]=sys.argv[1]+"F";
			infile = config['fileopts']['infolder']+'/'+sys.argv[1]
		if 'TAT' in config['prepro']['filter']:
			adapter='file:/home/labs/barkailab/LAB/scripts/felix_bcl2fastq/tat_adap.txt'		
			#cmd="/home/labs/barkailab/LAB/scripts/felix_bcl2fastq/cutadapt -b %s -B %s -j 16 -o %sF_R1_001.fastq.gz -p %sF_R2_001.fastq.gz %s_R1_001.fastq.gz %s_R2_001.fastq.gz -O 10 --pair-filter=any --max-n 0.8 --action=mask"
			cmd="cutadapt -b %s -B %s -j 16 -o %sF_R1_001.fastq.gz -p %sF_R2_001.fastq.gz %s_R1_001.fastq.gz %s_R2_001.fastq.gz -O 10 --pair-filter=any -m 15 --action=trim" # m was 18
			print(cmd %(adapter,adapter,infile,infile,infile,infile))
			psout=subprocess.check_output(cmd %(adapter,adapter,infile,infile,infile,infile),shell = True, stderr=subprocess.STDOUT) #
			print(psout)
			if (type(psout)==bytes):
				psout=psout.decode('utf-8');	
			adap_stats=re.findall('(?<=Trimmed\: ).*?(?= times)|(?<=read pairs processed\:)\s*[0-9,]*|(?<=passing filters\)\:).*(?= \()',psout)
			print(adap_stats)
			adap_line='\t'.join(['%.1f']*len(adap_stats))
			adap_stats = [adap_stat.replace(',','') for adap_stat in adap_stats]		
			adap_stats=map(float,adap_stats)
			adap_stats=numpy.asarray(list(adap_stats))
			adap_stats[1:]=100*adap_stats[1:]/adap_stats[0]
			os.system('echo %s >> %s/tat_stats.txt'% (sys.argv[1]+'\t'+adap_line %tuple(adap_stats),config['fileopts']['tempfolder']))
			#psout=re.findall('(?<=too many N\:).*(?= \()',psout)
			psout=re.findall('(?<=too short\:).*(?= \()',psout)
			prepro=prepro+int(psout[0].replace(',','').replace(' ',''))
			sys.argv[1]=sys.argv[1]+"F";
			infile = config['fileopts']['infolder']+'/'+sys.argv[1]
		if 'ATAC' in config['prepro']['filter']:
			adapter='file:/home/labs/barkailab/LAB/scripts/felix_bcl2fastq/atac_adap.txt'		
			#cmd="/home/labs/barkailab/LAB/scripts/felix_bcl2fastq/cutadapt -b %s -B %s -j 16 -o %sF_R1_001.fastq.gz -p %sF_R2_001.fastq.gz %s_R1_001.fastq.gz %s_R2_001.fastq.gz -O 10 --pair-filter=any --max-n 0.8 --action=mask"
			cmd="cutadapt -b %s -B %s -j 16 -o %sF_R1_001.fastq.gz -p %sF_R2_001.fastq.gz %s_R1_001.fastq.gz %s_R2_001.fastq.gz -O 10 --pair-filter=any -m 15 --action=trim" # m was 18
			print(cmd %(adapter,adapter,infile,infile,infile,infile))
			psout=subprocess.check_output(cmd %(adapter,adapter,infile,infile,infile,infile),shell = True, stderr=subprocess.STDOUT) #
			print(psout)
			if (type(psout)==bytes):
				psout=psout.decode('utf-8');	
			adap_stats=re.findall('(?<=Trimmed\: ).*?(?= times)|(?<=read pairs processed\:)\s*[0-9,]*|(?<=passing filters\)\:).*(?= \()',psout)
			print(adap_stats)
			adap_line='\t'.join(['%.1f']*len(adap_stats))
			adap_stats = [adap_stat.replace(',','') for adap_stat in adap_stats]		
			adap_stats=map(float,adap_stats)
			adap_stats=numpy.asarray(list(adap_stats))
			adap_stats[1:]=100*adap_stats[1:]/adap_stats[0]
			os.system('echo %s >> %s/tat_stats.txt'% (sys.argv[1]+'\t'+adap_line %tuple(adap_stats),config['fileopts']['tempfolder']))
			#psout=re.findall('(?<=too many N\:).*(?= \()',psout)
			psout=re.findall('(?<=too short\:).*(?= \()',psout)
			prepro=prepro+int(psout[0].replace(',','').replace(' ',''))
			sys.argv[1]=sys.argv[1]+"F";
			infile = config['fileopts']['infolder']+'/'+sys.argv[1]
	outfile = config['fileopts']['outfolder']+'/'+sys.argv[1]
	tempfile = config['fileopts']['tempfolder']+'/'+sys.argv[1]

	statfile = config['fileopts']['statfile']#'/home/labs/barkailab/felixj/scripts/dna/stats.xls'
	# two different bowtie commands depending on the file name of the read files - "_R1_001.fastq.gz" is the NovaSeq, _1 and _2 are downloaded from the internet.
	if 'R1' in config['fileopts']['ending']:
		cmd = '/home/labs/barkailab/felixj/bowtie2-2.3.5.1-source/bowtie2-2.3.5.1/bowtie2 %s -x %s -1 %s_R1_001.fastq.gz -2 %s_R2_001.fastq.gz | samtools sort -o %s.bam'
	else:
		cmd = 'bowtie2 %s -x %s -1 %s_1.fastq.gz -2 %s_2.fastq.gz | samtools sort -o %s.bam'
	print(cmd %(config['bowtie2']['args'],config['bowtie2']['indices'],infile,infile,tempfile))
	bt2out = subprocess.check_output(cmd %(config['bowtie2']['args'],config['bowtie2']['indices'],infile,infile,tempfile),shell = True, stderr=subprocess.STDOUT)
	
	# the pipeline has the ability to find the tagged ORF by loooking for DNA fragments that on the one end align to an ORF and on the other side to a tag (e.g. MNAse)
	if config.has_option('bowtie2','checkMNase'):
		cmd='/home/labs/barkailab/felixj/bowtie2-2.3.5.1-source/bowtie2-2.3.5.1/bowtie2 --quiet -p8 --very-sensitive --trim-to 30 --dovetail -x /home/labs/barkailab/felixj/RefData/yeast/Sc64KLac_cDna -1 %s_R1_001.fastq.gz -2 %s_R2_001.fastq.gz | samtools view -f1 -F2 -F4 -F8 | awk \'$3 ~ /%s/\' | cut -f7 | sort | uniq -c | sort -nr |head -1' 
		print(cmd % (infile,infile,config['bowtie2']['checkMNase']))
		checkOrf=subprocess.check_output(cmd % (infile,infile,config['bowtie2']['checkMNase']),shell = True, stderr=subprocess.STDOUT);
		checkOrf=checkOrf.decode('utf-8');
		checkOrf=checkOrf.replace('\n','')
	if (type(bt2out)==bytes):
		bt2out=bt2out.decode('utf-8');
	if config.has_section('prepro'):
		subprocess.check_output("rm %s_R*_001.fastq.gz" %(infile),shell = True, stderr=subprocess.STDOUT)
	
	#bowtie2 statistics are collected
	bt2out=list(map(int,re.findall('\d+(?= \(\d+)',bt2out)))
	tot_reads=bt2out[0]
	per_al=100*bt2out[2]/bt2out[0]
	mul_al=100*bt2out[3]/bt2out[0]
	dis_al=100*bt2out[4]/bt2out[0]
	mate_al=100*(bt2out[6]+bt2out[7])/bt2out[0]
	none_al=100*(bt2out[5]-(bt2out[6]+bt2out[7]))/(2*bt2out[0])
	al_pairs=bt2out[2]+bt2out[3]+bt2out[4];

	#a modified version of picard is used to mark PCR duplicates
	cmd_picard='java -jar /home/labs/barkailab/felixj/scripts/picard/build/libs/picard.jar MarkDuplicates %s I=%s.bam O=%s_wd.bam M=%s_dupl.tsv'
	print(cmd_picard%(config['picard']['args'],tempfile,tempfile,tempfile))
	picard_out=subprocess.check_output(cmd_picard%(config['picard']['args'],tempfile,tempfile,tempfile),shell = True, stderr=subprocess.STDOUT)

	if (type(picard_out)==bytes):
		picard_out=picard_out.decode('utf-8');
	
	#picard statisitcs are collected
	picard_out=re.search('(?<=PC duplicates )(\d+) size (\d+)\[',picard_out);
	if type(picard_out) is type(None):
		per_dup=0
		uni_pairs=al_pairs
		lib_size= 9999999
	else:
		per_dup=100*int(picard_out.group(1))/al_pairs
		uni_pairs=al_pairs-int(picard_out.group(1))
		lib_size=picard_out.group(2);
	
	#remove BAM files
	if os.path.isfile(tempfile+'_wd.bam'):
		os.system('rm '+tempfile+'.bam')
	else:
		os.system('mv '+tempfile+'.bam '+tempfile+'_wd.bam')
else: # this part of the pipeline is just run if bam files already exist so that the bowtie and picard statisitics get replaced by 0
	print('bam file exists')
	if config.has_option('bowtie2','checkMNase'):
		cmd='/home/labs/barkailab/felixj/bowtie2-2.3.5.1-source/bowtie2-2.3.5.1/bowtie2 --quiet -p8 --very-sensitive --trim-to 30 --dovetail -x /home/labs/barkailab/felixj/RefData/yeast/Sc64KLac_cDna -1 %s_R1_001.fastq.gz -2 %s_R2_001.fastq.gz | samtools view -f1 -F2 -F4 -F8 | awk \'$3 ~ /%s/\' | cut -f7 | sort | uniq -c | sort -nr |head -1' 
		print(cmd % (infile,infile,config['bowtie2']['checkMNase']))
		checkOrf=subprocess.check_output(cmd % (infile,infile,config['bowtie2']['checkMNase']),shell = True, stderr=subprocess.STDOUT);
		checkOrf=checkOrf.decode('utf-8');
		checkOrf=checkOrf.replace('\n','');		
	if config.has_section('prepro'):
		sys.argv[1]=sys.argv[1]+'F'
	tot_reads=0
	prepro=0
	per_al=0;
	mul_al=0
	dis_al=0
	mate_al=0
	none_al=0
	per_dup=0
	lib_size=0;
	al_pairs=0;
	outfile = config['fileopts']['outfolder']+'/'+sys.argv[1]
	tempfile = config['fileopts']['tempfolder']+'/'+sys.argv[1]
	statfile = config['fileopts']['statfile']
	
#median insertsize is calculated using samtools
cmd_insertsize='samtools stats -f 0x002 -F0x0400 %s_wd.bam | grep ^IS | cut -f2,3'
insertsize_out=subprocess.check_output(cmd_insertsize%(tempfile),shell = True, stderr=subprocess.STDOUT)
if (type(insertsize_out)==bytes):
	insertsize_out=insertsize_out.decode('utf-8');
insertsize_out=(re.findall('\d+(?=\n)',insertsize_out))
insertsize_out=list(map(int,insertsize_out))
if float(config['genomecov']['minLen'])<1 and float(config['genomecov']['minLen'])>0:
	print((sum(numpy.cumsum(insertsize_out)<(float(config['genomecov']['minLen'])*sum(insertsize_out)))))
insertsize_out=(sum(numpy.cumsum(insertsize_out)<(0.5*sum(insertsize_out)))) #calculate median size from histogram


## the pipeline has the option to run some matlab scripts to filter the aligned fragments using matlab... but its not really used (so just move to else)
if config.has_option('genomecov','preMat'):
	matCmd='matlab -nodisplay -r "filterBAM2 %s  ;exit"'%tempfile
	print(matCmd)
	matCmdOut=subprocess.check_output(matCmd,shell = True, stderr=subprocess.STDOUT)
	cmd = 'samtools view -h %s %s_wd.bam | awk -F "\\t" \'{if ((($9^2>=%s^2)&&($9^2<=%s^2)&&'+config['genomecov']['addconstraint']+') || ($1~"@HD|@PG") || (($1~"@SQ")&&($2~"NC_001")) )  print $_}\' | samtools view -c'
	woConstr_pairs=0;
else:
	## usually the pipeline align allr eads in between the defined sizes- however its possible to add additonal constrains based on the alignment scores ($_~"AS:i:0"||$_~"YS:i:0"))
	if config.has_option('genomecov','addconstraint'):
		cmd = 'samtools view -h %s %s_wd.bam | awk -F "\\t" \'{if ((($9^2>=%s^2)&&($9^2<=%s^2)) || ($1~"@"))  print $_}\' | grep -v "XS:"|samtools view -c';
		woConstr_pairs=subprocess.check_output(cmd %(config['genomecov']['sam_flag'],tempfile,config['genomecov']['minLen'],config['genomecov']['maxLen']),shell = True, stderr=subprocess.STDOUT)
		if (type(woConstr_pairs)==bytes):
			woConstr_pairs=woConstr_pairs.decode('utf-8');
		woConstr_pairs=int(woConstr_pairs)/2;
		#cmd = 'samtools view -h %s %s_wd.bam | awk -F "\\t" \'{if ((($9^2>=%s^2)&&($9^2<=%s^2)&&'+config['genomecov']['addconstraint']+') || ($1~"@HD|@PG") || (($1~"@SQ")&&($2~"NC_001")) )  print $_}\' |grep -v "XS:"| samtools view -c'
		cmd = 'samtools view -h %s %s_wd.bam | awk -F "\\t" \'{if ((($9^2>=%s^2)&&($9^2<=%s^2)&& ($_~"AS:i:0"||$_~"YS:i:0")) || ($1~"@HD|@PG") || (($1~"@SQ")))  print $_}\' | samtools view -c'
	else: ##this is the most commonly used command.
		cmd = 'samtools view -h %s %s_wd.bam | awk -F "\\t" \'{if ((($9^2>=%s^2)&&($9^2<=%s^2)) || ($1~"@"))  print $_}\' | samtools view -c'
## before the actual coverage calculation is done - first the number of reads that go into it is calcualted.

print(cmd %(config['genomecov']['sam_flag'],tempfile,config['genomecov']['minLen'],config['genomecov']['maxLen']))
uni_pairs=subprocess.check_output(cmd %(config['genomecov']['sam_flag'],tempfile,config['genomecov']['minLen'],config['genomecov']['maxLen']),shell = True, stderr=subprocess.STDOUT)
if (type(uni_pairs)==bytes):
	uni_pairs=uni_pairs.decode('utf-8');
uni_pairs=int(uni_pairs)/2


## now the same thing happens just with that instead of counting the reads in the end "samtools view -c" they are piped into genomeCOverage "| samtools view -b | %s %s -d  -ibam - | cut -f3 >%s.out'"
if config.has_option('genomecov','genomeCovBin'):
	genomeCovBin=config['genomecov']['genomeCovBin']
else:
	genomeCovBin='genomeCoverageBed'
if config.has_option('genomecov','preMat'):
		cmd = 'samtools view -h %s %s_wd.bam | awk -F "\\t" \'{if ((($9^2>=%s^2)&&($9^2<=%s^2)&&'+config['genomecov']['addconstraint']+') || ($1~"@HD|@PG") || (($1~"@SQ")&&($2~"NC_001")))  print $_}\' |samtools view -b | %s %s -d  -ibam - | cut -f3 >%s.out'
else:
	if config.has_option('genomecov','addconstraint'):
		#cmd = 'samtools view -h %s %s_wd.bam | awk -F "\\t" \'{if ((($9^2>=%s^2)&&($9^2<=%s^2)&&'+config['genomecov']['addconstraint']+') || ($1~"@HD|@PG") || (($1~"@SQ")&&($2~"NC_001")))  print $_}\' | grep -v "XS:"|samtools view -b | %s %s -d  -ibam - | cut -f3 >%s.out'
		cmd = 'samtools view -h %s %s_wd.bam | awk -F "\\t" \'{if ((($9^2>=%s^2)&&($9^2<=%s^2)&& ($0~"AS:i:0"||$0~"YS:i:0")) || ($1~"@HD|@PG") || (($1~"@SQ")))  print $_}\' |samtools view -b | %s %s -d  -ibam - | cut -f3 >%s.out'
	else:
		cmd = 'samtools view -h %s %s_wd.bam | awk -F "\\t" \'{if ((($9^2>=%s^2)&&($9^2<=%s^2)) || ($1~"@"))  print $_}\' | samtools view -b | %s %s -d  -ibam - | cut -f3 >%s.out'
print(cmd %(config['genomecov']['sam_flag'],tempfile,config['genomecov']['minLen'],config['genomecov']['maxLen'],genomeCovBin,config['genomecov']['args'],outfile))
gcout=subprocess.check_output(cmd %(config['genomecov']['sam_flag'],tempfile,config['genomecov']['minLen'],config['genomecov']['maxLen'],genomeCovBin,config['genomecov']['args'],outfile),shell = True, stderr=subprocess.STDOUT)

## all the collected stats are written into a stats file that combines all samples.
if config.has_section('prepro'):
	statLine='%s\t%s\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%s\t%.1f\t%.1f\t%s\t%s';
	statline = statLine%(sys.argv[1],prepro,tot_reads,per_al,mul_al,dis_al,mate_al,none_al,al_pairs,per_dup,uni_pairs,lib_size,insertsize_out)
else:
	statLine='%s\t%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%s\t%.1f\t%.1f\t%s\t%s';
	statline = statLine%(sys.argv[1],tot_reads,per_al,mul_al,dis_al,mate_al,none_al,al_pairs,per_dup,uni_pairs,lib_size,insertsize_out)
if config.has_option('bowtie2','checkMNase'):
	statline=statline+'\t'+checkOrf;
if config.has_option('genomecov','addconstraint'):
	print(type(woConstr_pairs))
	print(woConstr_pairs)
	statline=statline+'\t%.0f'%(woConstr_pairs);
subprocess.call('echo "'+statline+'">>'+statfile, shell =True)


# df = pandas.read_excel(stat_file)
# dfnew=df.append(pandas.Series([sys.argv[1],bt2out[0],bt2out[2],bt2out[3],bt2out[4],bt2out[6]+bt2out[7],bt2out[5]-(bt2out[6]+bt2out[7]),n_dupl,insertsize_out],index=df.columns ),ignore_index=True)
# dfnew.to_excel(stat_file, engine='xlsxwriter')
