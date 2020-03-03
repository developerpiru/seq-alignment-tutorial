# Cloud servers
These are a collection of scripts and instructions to leverage cloud computing platforms for NGS analysis (like RNA-seq)

## Table of contents

+ [Cloud platforms](https://github.com/developerpiru/cloudservers#cloud-platforms)

  + [Google Cloud Platform (GCP)](https://github.com/developerpiru/cloudservers#google-cloud-platform)
+ [Getting started with GCP](https://github.com/developerpiru/cloudservers#getting-started-with-gcp)

## Cloud platforms

### Google Cloud Platform (GCP)

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
  
## RNA-seq analysis with your VM

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
Download the full assembly ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gzor only download the region of interest depending on your experiment. 

  ```
  wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  gunzip -k Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  ```
  This may take a while.

Second, you need a gtf which holds gene structure information. You can access release 88 for the human genome here: http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/ENSEMBL/homo_sapiens/ENSEMBL.homo_sapiens.release-83/Homo_sapiens.GRCh38.83.gtf

Now make a new folder in your home directory to store the generated reference genome:
  ```
  cd $HOME
  mkdir STARgenome
  ```

Now you can enter the command to generate the genome:
  ```
  STAR --runThreadN 64 --runMode genomeGenerate --genomeDir $HOME/STARgenome --genomeFastaFiles $HOME/STARgenomefiles/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile makegenomefiles/Homo_sapiens.GRCh38.83.gtf --sjdbOverhang 100
  ```
  Where ```--runThreadN 64``` specificies the number of CPU threads you want to use (general rule is CPU cores x2 or CPU cores x4)
  
This will take a while to generate.



