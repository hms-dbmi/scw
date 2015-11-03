# HSCI/Catalyst Single-cell RNA-Seq Workshop, November 4th, 2015

This repository contains R markdown files corresponding to the tutorial exercises that we will be running on Wednesda afternoon.

The code and output can be viewed at the accompanying website, http://hms-dbmi.github.io/scw

To begin, login to your training accounts on Orchestra

```{bash}
ssh -X trainingNN@orchestra.med.harvard
```
where `NN` should be replaced by your unique user number, which will be provided to you on the day.

Then start up an interactive session with a bash shell using the following command

```{bash}
bsub -n 2 -Is -q interactive bash
```

Now to download the material in the repository
```{bash}
git clone https://github.com/hms-dbmi/scw.git
cd scw/scw2015
```

Next we will run a small setup script to create the correct environment variables
```{bash}
source setup.sh
```

Now you can change directory into the subdirectory corresponding to the first exercise (alignment), and start `R` by executing the command
```{bash}
R
```

