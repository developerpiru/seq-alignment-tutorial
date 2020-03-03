# Cloud servers
These are a collection of scripts and instructions to leverage cloud computing platforms for NGS analysis (like RNA-seq)

## Table of contents

+ [Cloud platforms](https://github.com/developerpiru/BEAVR#cloud-platforms)

  + [Google Cloud Platform (GCP)](https://github.com/developerpiru/BEAVR#google-cloud-platform)
  + [Amazon Web Services (AWS)](https://github.com/developerpiru/BEAVR#amazon-web-services)
  + [Microsoft Azure](https://github.com/developerpiru/BEAVR#microsoft-azure)
+ [Getting started with GCP](https://github.com/developerpiru/BEAVR#getting-started-with-gcp)

## Cloud platforms

### Google Cloud Platform (GCP)

### Amazon Web Services (AWS)

### Microsoft Azure

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

4.
