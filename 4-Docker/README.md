# Introduction to Building Docker for BERLIN

### **Overview:**

In the dynamic realm of single-cell RNA sequencing (scRNA-seq), the analysis of intricate biological data holds paramount significance for unveiling the nuances of cellular heterogeneity and genetic expression. To simplify this intricate process and streamline scRNAseq data analysis, we present the Building Docker for BERLIN project.

### **BERLIN: Basic Exploration of single-cell RNAseq data and LINeages**

<https://github.com/shawlab-moffitt/BERLIN>

BERLIN, an acronym for "Basic Exploration of single-cell RNAseq data and LINeages," offers a comprehensive workflow for effective scRNAseq data analysis. Leveraging R, the Seurat package, and essential tools like SeuratDisk, pathwork, dplyr, Single R, and Celldex, this workflow is encapsulated within a Docker container, enhancing the deployment and reproducibility of scRNAseq analysis.

### The benefits of Docker for BERLIN are threefold:

1.  **Isolation and Reproducibility:** Docker containers encompass the entire environment, including dependencies, libraries, and tools, ensuring analysis consistency across diverse systems and environments. This isolation guarantees the reproducibility of results, making it easier for researchers to share their work and for others to replicate it.

2.  **Portability**: Docker containers are highly portable, capable of running on any system supporting Docker. This portability facilitates the effortless sharing of analysis workflows among collaborators, irrespective of their underlying operating systems or infrastructure, thereby streamlining collaboration.

3.  **Efficient Environment Setup**: Docker simplifies the setup of complex analysis environments, alleviating researchers from the time-consuming and error-prone task of manual software installation and configuration. With Docker, the entire environment is defined in code, enabling consistent and hassle-free setup.

### **Running the Docker Container Using Docker Desktop**:

To embark on Building Docker for BERLIN, we encourage you to install Docker Desktop on your system. Once you've installed Docker Desktop, you can build the Docker image and run the container with the following simple steps:

#### **Single Cell RNAseq Pipeline Docker Setup**:

**1. Install Docker Desktop**: Download and install Docker Desktop for your operating system from the official Docker website ([https://www.docker.com/products/docker-desktop).](https://www.docker.com/products/docker-desktop).)

**2. Build the Docker Image**:

-   Clone the Docker project repository from its source.

-   Open Windows command prompt and navigate to the directory containing the Dockerfie.

```         
cd C:\Users\Administrator\Desktop\BERLIN_02\2-Single_Cell_RNAseq_Pipeline_docker\Single_Cell_RNAseq_Pipeline_docker
```

-   Use the `docker build` command to build the Docker image. For example:

```         
   docker build -t single_cell_rnaseq_analysis_app .
```

**3. Run the Docker Container**:

-   Once the Docker image is successfully built, you can run it using Docker Desktop. For example, click on the **'run'** **icon** to start the container, and **all output files** will be generated inside the **'app' folder** of the downloaded repository.

    ![](https://github.com/chingyaousf/Introduction-to-Building-Docker-for-BERLIN-Pipeline/blob/main/data/Docker%20Desktop%20images_04.png?raw=true)

    ![](https://github.com/chingyaousf/Introduction-to-Building-Docker-for-BERLIN-Pipeline/blob/main/data/Docker%20Desktop%20images_05.png?raw=true)

#### **Shiny Visualization Applications Docker Setup (3 Containers, One Example Included):**

**1. Install Docker Desktop**: Download and install Docker Desktop for your operating system from the official Docker website ([https://www.docker.com/products/docker-desktop).](https://www.docker.com/products/docker-desktop).)

**2. Build the Docker Image**:

-   Clone the Docker project repository from its source.

-   Open Windows command prompt and navigate to the directory containing the Dockerfie.

```         
cd C:\Users\Administrator\Desktop\BERLIN_02\3-R_Shiny_Viz_Applications_docker\DRPPM_EASY_Shiny_App_docker
```

-   Use the `docker build` command to build the Docker image. For example:

```         
   docker build -t drppm_easy_shiny_app .
```

**3. Run the Docker Container**:

-   Once the Docker image is successfully built, you can run the image using Docker Desktop. For instance, click 3838:3838 to display web page as following.

    ![](https://github.com/chingyaousf/Introduction-to-Building-Docker-for-BERLIN/blob/main/data/Docker%20Desktop%20images.jpg?raw=true)

![](https://github.com/chingyaousf/Introduction-to-Building-Docker-for-BERLIN/blob/main/data/Docker%20Desktop%20images_02.jpg?raw=true)

![](https://github.com/chingyaousf/Introduction-to-Building-Docker-for-BERLIN/blob/main/data/Docker%20Desktop%20images_03.jpg?raw=true)

**4. Use Docker Compose for Building and Running**:

-   Create a Docker Compose YAML file in the Docker project directory (already included in the Docker file).

-   Inside the YAML file, specify the services, image names, ports, and other configurations.

-   Use **`docker-compose build`** to build the image based on the YAML file. You can run this command on Windows command prompt whenever the Dockerfile is modified to update the image.

-   Use **`docker-compose up`** to run the container based on the built image. This is also used to start the container whenever needed.

-   To remove the container and associated resources, use **`docker-compose down`**.

The Docker container is configured to run the R and Shiny app, ensuring that you can efficiently explore the intricacies of single-cell RNA sequencing with ease and consistency. This setup allows for convenient access to the R and Shiny app, which output can be accessed on your local system inside app folder or via the specified port. Remember that you can run `docker-compose build` or `docker-compose up` every time you modify the Dockerfile to update the image and keep your analysis environment up to date. Please refer to the **Docker Desktop Learning Center** and **Docker Compose documentation** for comprehensive guidance on using these tools effectively.

## Reference:

### **BERLIN:**

<https://github.com/shawlab-moffitt/BERLIN>

## Blog:

<https://ssidmarine.wordpress.com/2023/11/01/introduction-to-building-docker-for-berlin/>
