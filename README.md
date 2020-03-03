# Cloud servers
These are a collection of scripts and instructions to leverage cloud computing platforms for NGS analysis (like RNA-seq)

## Table of contents

+ [Cloud platforms](https://github.com/developerpiru/cloudservers#cloud-platforms)

  + [Google Cloud Platform (GCP)](https://github.com/developerpiru/cloudservers#google-cloud-platform-gcp)
    + [Getting started with GCP](https://github.com/developerpiru/cloudservers#getting-started-with-gcp)
    + [Creating a new compute instance](https://github.com/developerpiru/cloudservers#creating-a-new-compute-instance)
    + [Creating a storage bucket](https://github.com/developerpiru/cloudservers#creating-a-storage-bucket)
    + How to setup a newly created VM](https://github.com/developerpiru/cloudservers#how-to-setup-a-newly-created-vm)
      + [Setting up GCSFuse](https://github.com/developerpiru/cloudservers#setting-up-gcsfuse)
      + [Install Miniconda](https://github.com/developerpiru/cloudservers#install-miniconda)
      + [Install Basemount](https://github.com/developerpiru/cloudservers#installing-basemount)
    + RNA-seq analysis with your VM](https://github.com/developerpiru/cloudservers#rna-seq-analysis-with-your-vm)     
      + [Setting up STAR aligner](https://github.com/developerpiru/cloudservers#setting-up-star-aligner)
      + [Preparing your read files](https://github.com/developerpiru/cloudservers#preparing-your-read-files)
      + [Get read counts with STAR](https://github.com/developerpiru/cloudservers#get-read-counts-with-star)

---

## Cloud platforms

### Google Cloud Platform (GCP)

Create a Google account if you don't already have one and then go to https://console.cloud.google.com/ to setup your cloud account. 

## Getting started with GCP

### Creating a new compute instance

1. From the [GCP Dashboard](https://console.cloud.google.com/home/dashboard), click on the menu button (1), then hover over **Computer Engine** to display the pop-up menu (2), then click on **VM instances** (3)
  ![Image of GCP Dashboard](Screenshots/GCP/createinstance.png)
  
2. This is where all of your existing VMs will be listed. To create a new VM, click **Create Instance** at the top.
  ![Image of GCP Dashboard](Screenshots/GCP/createinstance-2.png)

3. The page will ask you to specify details for your new VM instance.
  ![Image of GCP Dashboard](Screenshots/GCP/createinstance-3.png)
  
  You need to specify these details:
  + Instance name
  + Region and zone: Select **us-east1** for Region because it is closest to you and the cheapest cost.
  + Machine configuration: select a prebuilt configuration or customize your own CPU and memory allotment. 
  + Click the **Change** button under **Boot disk** to specify the hard drive and OS configuration
![Image of GCP Dashboard](Screenshots/GCP/createinstance-4.png)
    + Select Ubuntu for the **Operating system**
    + Select Ubuntu 18.04 LTS for **Version**
    + Set the size for 500 GB
    
  + Under **Identity and API access**, select **Allow full access to all Cloud APIs**
  + Under **Firewall**, select **Allow HTTP traffic** and **Allow HTTPS traffic**
  + Click the **Create** button at the button when you are done.

4. This will return you to the **VM instances** page where you will see your newly created instance. The green icon indicates the VM is currently running.
![Image of GCP Dashboard](Screenshots/GCP/createinstance-5.png)

5. To connect to the VM, click on the **SSH** button (after ensuring it is on). You can use the menu button to Start, Stop, Reset, or Delete a VM.

### Creating a storage bucket

Think of a storage bucket as a cloud storage drive where you can save your files. You can use it as a shared folder to transfer files to and from your local computer to a VM. 

1. From the [GCP Dashboard](https://console.cloud.google.com/home/dashboard), click on the menu button, then click on **Storage**. 

2. This will open the storage browser and list any buckets you have (similar to the VM Instances page. Click **Create Bucket** at the top.

3. This will open a new page where you must enter details for your new bucket.
![Image of GCP Dashboard](Screenshots/GCP/createbucket-1.png)
  + You must specify a globally unique name (meaning the name must not be taken by anyone else)
  + For **Location type**, select Region
  + For **Location**, select the **same region as the VM instances you have created!**
  + Leave all other settings as defaults and click **Create** at the bottom.

4. Your newly created bucket will be listed. Click on its name to view its contents (which will be empty right now).

5. You can upload files by using the buttons provided or by dragging files into the browser from your computer.

### How to setup a newly created VM

#### Setting up GCSFuse
A newly created VM only has a basic install of Ubuntu and you must now install any software you need. 

The first step is to install [gcsfuse](https://github.com/GoogleCloudPlatform/gcsfuse) on your new VM. Gcsfuse allows you to connect your VM to a storage bucket you own. This will "mount" your storage bucket as a virtual drive in your VM so you can transfer files to and from it. 

Copy and paste these commands to setup gcsfuse:
  ```
  export GCSFUSE_REPO=gcsfuse-`lsb_release -c -s`
  echo "deb http://packages.cloud.google.com/apt $GCSFUSE_REPO main" | sudo tee /etc/apt/sources.list.d/gcsfuse.list
  curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo apt-key add -
  sudo apt-get update
  sudo apt-get install gcsfuse
  mkdir mountfolder
  ```

You now have gcsfuse installed and have created a folder called "mountfolder". This is where you will mount your storage bucket.

Enter this command to mount your storage bucket:
  ```
  gcsfuse myBucketName mountfolder
  ```
  Where ```myBucketName``` is the name of the storage bucket you created.

#### Install Miniconda

Miniconda is required to download and install some bioinformatics tools. You can install it using these commands:
  ```
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh
  ```
  Then follow the prompts to complete installation.
  
You may have to open a new terminal window in order for the ```conda``` command to be recognized.

Try entering this command:
  ```
  conda --version
  ```
  You should see the version number of conda if it installed correctly.
  
#### Installing Basemount

Like gcsfuse, basemount has their own software to let you interface with your basepsace account. This will allow you to mount your basespace account as a virtual drive and transfer your read files.

Enter these commands to install basemount:
  ```
  sudo bash -c "$(curl -L https://basemount.basespace.illumina.com/install)"
  mkdir basemountfolder
  ```
  
You now have basemount installed and have created a folder called "basemountfolder". This is where you will mount your basemount account.

Enter this command to mount your storage bucket:
  ```
  basemount basemountfolder
  ```
  You will now be prompted to login to your basespace account and once you do, you can access your files.

---

## RNA-seq analysis with your VM

### Setting up STAR aligner

Once you have the basics setup, you can begin your RNA-seq analysis. 

First you need to install STAR aligner. 
  ```
  conda install STAR
  ```

Then you need to download some required reference files in order to build a reference genome for your organism of interest (e.g. human, mouse)

Make a new directory to store these files and move into that folder:
  ```
  cd $HOME
  mkdir STARgenomefiles
  cd STARgenomefiles
  ```

First, you need a "primary assembly file" in FASTA format. You can access it here: ftp://ftp.ensembl.org/pub/ 
Look for the latest release, e.g. for humans: ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/
Download the full assembly ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz or only download the region of interest depending on your experiment. 

  ```
  wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  gunzip -k Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  ```
  This may take a while.

Second, you need a gtf which holds gene structure information. Make sure the release version is the same as your primary assembly file above. You can access release 99 for the human genome here: ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz

  ```
  wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
  gunzip -k homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
  ```

Now make a new folder in your home directory to store the generated reference genome:
  ```
  cd $HOME
  mkdir STARgenome
  ```

Now you can enter the command to generate the genome:
  ```
  STAR --runThreadN 64 --runMode genomeGenerate --genomeDir $HOME/STARgenome --genomeFastaFiles $HOME/STARgenomefiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile makegenomefiles/Homo_sapiens.GRCh38.99.gtf --sjdbOverhang 100
  ```
  Where ```--runThreadN 64``` specificies the number of CPU threads you want to use (general rule is CPU cores x2 or CPU cores x4)
  
This will take a while to generate.

### Preparing your read files

Mount your basemount folder and download your files to your VM:
  ```
  basemount basemountfolder
  ```

Create a new folder to store your downloaded fastq files
  ```
  mkdir myfastqfiles
  ```

Once you have downloaded your files, extract them:
  ```
  gunzip -k *.gz
  ```

Then you need to concatenate files run in different names for the same sample:
You can use wildcards like this:
  ```
  cat WTA-1*.fastq > WTA-1.fastq
  cat WTA-2*.fastq > WTA-2.fastq
  cat WTA-3*.fastq > WTA-3.fastq
  
  cat KO4A-1*.fastq > KO4A-1.fastq
  cat KO4A-2*.fastq > KO4A-2.fastq
  cat KO4A-3*.fastq > KO4A-3.fastq
  ```
  
### Get read counts with STAR

You can run this command to start alignments and generating read counts with STAR. Modify these parameters for your experiment:
  + This example shows how to do this for 6 samples with fastq files named: ```WTA-1.fastq```, ```WTA-2.fastq```, ```WTA-3.fastq```, ```KOA-1.fastq```, ```KOA-2.fastq```, ```KOA-3.fastq```. 
  + The fastq files are stored in ```~/mount/RNAseq-aug2018/``` which is specified with ```--readFilesIn```.
  + The location of the genome you generated is specified by ```--genomeDir```
  + The location of the gtf file you downloaded is specififed by ```--sjdbGTFfile```
  + The output location and filename prefix is specified by ```--outFileNamePrefix```

```
#WTA-1.fastq:
STAR --runThreadN 96 --genomeDir ~/mount/STARgenome --readFilesIn ~/mount/RNAseq-aug2018/WTA-1.fastq --sjdbGTFfile ~/mount/STARfiles/makegenomefiles/Homo_sapiens.GRCh38.92.gtf --sjdbOverhang 100 --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix ~/mount/RNAseq-aug2018/20180814-readcounts/WTA-1 --genomeLoad LoadAndKeep

#WTA-2.fastq:
STAR --runThreadN 96 --genomeDir ~/mount/STARgenome --readFilesIn ~/mount/RNAseq-aug2018/WTA-2.fastq --sjdbGTFfile ~/mount/STARfiles/makegenomefiles/Homo_sapiens.GRCh38.92.gtf --sjdbOverhang 100 --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix ~/mount/RNAseq-aug2018/20180814-readcounts/WTA-2

#WTA-3.fastq:
STAR --runThreadN 96 --genomeDir ~/mount/STARgenome --readFilesIn ~/mount/RNAseq-aug2018/WTA-3.fastq --sjdbGTFfile ~/mount/STARfiles/makegenomefiles/Homo_sapiens.GRCh38.92.gtf --sjdbOverhang 100 --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix ~/mount/RNAseq-aug2018/20180814-readcounts/WTA-3

#KO4A-1.fastq:
STAR --runThreadN 96 --genomeDir ~/mount/STARgenome --readFilesIn ~/mount/RNAseq-aug2018/KO4A-1.fastq --sjdbGTFfile ~/mount/STARfiles/makegenomefiles/Homo_sapiens.GRCh38.92.gtf --sjdbOverhang 100 --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix ~/mount/RNAseq-aug2018/20180814-readcounts/KO4A-1

#KO4A-2.fastq:
STAR --runThreadN 96 --genomeDir ~/mount/STARgenome --readFilesIn ~/mount/RNAseq-aug2018/KO4A-2.fastq --sjdbGTFfile ~/mount/STARfiles/makegenomefiles/Homo_sapiens.GRCh38.92.gtf --sjdbOverhang 100 --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix ~/mount/RNAseq-aug2018/20180814-readcounts/KO4A-2

#KO4A-3.fastq:
STAR --runThreadN 96 --genomeDir ~/mount/STARgenome --readFilesIn ~/mount/RNAseq-aug2018/KO4A-3.fastq --sjdbGTFfile ~/mount/STARfiles/makegenomefiles/Homo_sapiens.GRCh38.92.gtf --sjdbOverhang 100 --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix ~/mount/RNAseq-aug2018/20180814-readcounts/KO4A-3
```

