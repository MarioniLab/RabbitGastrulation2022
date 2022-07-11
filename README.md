# RabbitGastrulation2022

This repository contains the functions and scripts used as part of the Ton et al. 2022 paper on rabbit gastrulation. 



> **Rabbit Development as a Model for Single Cell Comparative Genomics**
>
> Mai-Linh Ton<sup>1,2,\*</sup> Daniel Keitley<sup>3,*</sup>, Bart Theeuwes<sup>2</sup>, Carolina Guibentif<sup>4</sup>, Jonas Ahnfelt-Rønne<sup>5</sup>, Thomas Kjærgaard Andreassen<sup>5</sup>, Fernando J. Calero-Nieto<sup>2</sup>, Ivan Imaz-Rosshandler<sup>6</sup>, Blanca Pijuan-Sala<sup>7</sup>, Jennifer Nichols<sup>2</sup>, Èlia Benito-Gutiérrez<sup>3</sup>, John Marioni<sup>7,8,9</sup>, Berthold Göttgens<sup>1,2</sup>



This project relies on both R and python code packaged into [scrabbitr](https://github.com/dkeitley/scrabbitr) and [scrabbitpy](https://github.com/dkeitley/scrabbitpy) respectively. We recommend installing these packages within R and python conda environments. The environments originally used in our analysis can be found in the [envs](envs/) directory.  

##### Installing scrabbitr

```R
# Inside R session
devtools::install_github("https://github.com/dkeitley/scrabbitr")
```

##### Installing scrabbit

```bash
git clone https://github.com/dkeitley/scrabbitpy.git
cd scrabbitpy
pip install .
```



