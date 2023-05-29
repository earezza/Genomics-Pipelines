#!/bin/bash
# After logging in, enter -> cd projects/def-jdilwort/$USER/
# Then the script below will prep everything and may take several minutes
# Note that only mm10 and hg38 genome assemblies are downloaded automatically, others require manual downloading

date

# --- Download GitHub ---
#git clone https://github.com/Dilworth-Laboratory/Genomics-Pipelines.git
cp Genomics-Pipelines/Process-NGS-RawReads/ngs_processing_pipeline.py .
cp Genomics-Pipelines/Process-NGS-RawReads/ngs_pipeline_*.sh .


# --- Download external software ---
wget -O picard.jar https://github.com/broadinstitute/picard/releases/download/3.0.0/picard.jar
git clone https://github.com/FredHutch/SEACR.git
wget -O gopeaks https://github.com/maxsonBraunLab/gopeaks/releases/download/v1.0.0/gopeaks-linux-amd64
chmod +x gopeaks


# === Create virtualenvironment ===
if test -d ngsENV/; then 
	echo "virtualenv ngsENV already exists..."; 
else 
	echo "Creating virtualenv" 
	module load python/3.9
	virtualenv --no-download ngsENV
	source ngsENV/bin/activate
	pip install --no-index --upgrade pip
	pip install -r Genomics-Pipelines/Process-NGS-RawReads/cc_requirements.txt
fi


# === Download reference genomes ===
if test -d ../Reference_Files/; then
	echo "Shared ../Reference_Files/ folder found..."; 
	
	# --- Look for bowtie2 files and download if not found ---
	if test -d ../Reference_Files/bowtie2/; then
		
		echo "bowtie2/ folder found...";
		
		# Look for mm10 and download if not found
		if test -f ../Reference_Files/bowtie2/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome.fa; then 
			echo "mm10 found"; 
		else 
			cd ../Reference_Files/bowtie2/
			echo "Downloading mm10";
			wget -O Mus_musculus_UCSC_mm10.tar.gz http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
			tar -xzvf Mus_musculus_UCSC_mm10.tar.gz 
			cd ../../$USER/
		fi
		
		# Look for hg38 and download if not found
		if test -f ../Reference_Files/bowtie2/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome.fa; then 
			echo "hg38 found"; 
		else 
			cd ../Reference_Files/bowtie2/
			echo "Downloading hg38"
			wget -O Homo_sapiens_UCSC_hg38.tar.gz http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
			tar -xzvf Homo_sapiens_UCSC_hg38.tar.gz 
			cd ../../$USER/
		fi
	else
		echo "Downloading bowtie2 genome index files...";
		mkdir ../Reference_Files/bowtie2/
		chmod +777 ../Reference_Files/bowtie2/
		cd ../Reference_Files/bowtie2/
		
		echo "Downloading mm10";
		wget -O Mus_musculus_UCSC_mm10.tar.gz http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
		tar -xzvf Mus_musculus_UCSC_mm10.tar.gz
		
		echo "Downloading hg38"
		wget -O Homo_sapiens_UCSC_hg38.tar.gz http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
		tar -xzvf Homo_sapiens_UCSC_hg38.tar.gz
		cd ../../$USER/
	fi
	
	# --- Look for hisat2 files and download if not found ---
	if test -d ../Reference_Files/hisat2/; then
		
		echo "hisat2/ folder found...";
		
		# Look for mm10 and download if not found
		if test -f ../Reference_Files/hisat2/mm10/genome.1.ht2; then 
			echo "mm10 found"; 
		else 
			cd ../Reference_Files/hisat2/
			echo "Downloading mm10";
			wget -O mm10_genome.tar.gz https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz
			tar -xzvf mm10_genome.tar.gz 
			cd ../../$USER/
		fi
		
		# Look for hg38 and download if not found
		if test -f ../Reference_Files/hisat2/hg38/genome.1.ht2; then 
			echo "hg38 found"; 
		else 
			cd ../Reference_Files/bowtie2/
			echo "Downloading hg38"
			wget -O hg38_genome.tar.gz https://genome-idx.s3.amazonaws.com/hisat/hg38_genome.tar.gz
			tar -xzvf hg38_genome.tar.gz 
			cd ../../$USER/
		fi
	else
		echo "Downloading hisat2 genome index files...";
		mkdir ../Reference_Files/hisat2/
		chmod +777 ../Reference_Files/hisat2/
		cd ../Reference_Files/hisat2/
		
		echo "Downloading mm10";
		wget -O mm10_genome.tar.gz https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz
		tar -xzvf mm10_genome.tar.gz
		
		echo "Downloading hg38"
		wget -O hg38_genome.tar.gz https://genome-idx.s3.amazonaws.com/hisat/hg38_genome.tar.gz
		tar -xzvf hg38_genome.tar.gz
		
		cd ../../$USER/
	fi
	
else 
	echo "Shared Reference_Files/ not found...creating now...";
	mkdir ../Reference_Files/
	chmod +777 ../Reference_Files/
	
	echo "Downloading bowtie2 genome index files...";
	mkdir ../Reference_Files/bowtie2/
	chmod +777 ../Reference_Files/bowtie2/
	cd ../Reference_Files/bowtie2/
	
	echo "Downloading mm10";
	wget -O Mus_musculus_UCSC_mm10.tar.gz http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
	tar -xzvf Mus_musculus_UCSC_mm10.tar.gz
	
	echo "Downloading hg38"
	wget -O Homo_sapiens_UCSC_hg38.tar.gz http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
	tar -xzvf Homo_sapiens_UCSC_hg38.tar.gz
	
	cd ../
	
	echo "Downloading hisat2 genome index files...";
	mkdir hisat2/
	chmod +777 hisat2/
	cd hisat2/
	
	echo "Downloading mm10";
	wget -O mm10_genome.tar.gz https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz
	tar -xzvf mm10_genome.tar.gz
	
	echo "Downloading hg38"
	wget -O hg38_genome.tar.gz https://genome-idx.s3.amazonaws.com/hisat/hg38_genome.tar.gz
	tar -xzvf hg38_genome.tar.gz

	cd ../../$USER/
fi
chmod +777 -R ../Reference_Files/*

date
