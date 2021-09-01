# Tissue/Gene specificity Scoring

## Introduction
In this project, the goal is to find the most tissue specific genes.
In order to do so, we use several scores that will allow us to rate that specificity.
There are mainly 2 types of scores: the generic scores that hints you whether a gene can have high tissue specificity by rating the gene with one value, and the specific scores that give you one score per tissue per gene.

We have used RNA-seq data in TPM format (or pTPM), note that counts data can be used too.

## Pre - Installation steps:
- Make sure you have python
- Make sure you have pip package
- Make sure you clone this package

## Install requirements:
```bash
pip install -r requirements.txt
```

## What you need to run:

1. To get generic scores, run the file `run-generic.py`
2. To get specific scores, first run the file `run-specs.py` if you want to include the SPECs score in your analysis (it's very computationally demanding, I performed it on a computer with 538 GB of ram and 128 cores), then run the file `run-specific.py`.

You will get outputs in the dataframes and the distributions folders.

## Workflow of the code

## Part 1. Filtering post computation of scores
### Step 1: Load and format the data
### Step 2: Compute scores (generic and specific)
### Step 3: Output distribution of data
### Step 4: Filtering
- Filter on subcellular location (genes) for generic and specific scores (see list)
- Filter on tissues for specific scores (highly specific ones)
### Step 5: Heatmaps & Intersections

## Part 2. Filtering pre- and post-computation of scores
### Step 1: Load and format the data
### Step 2: Filtering
- Filter on tissues for generic and specific scores (see list)
### Step 3: Compute scores (generic and specific)
### Step 4: Output distribution of data
### Step 5: Filtering
- Filter on subcellular location (genes) for generic and specific scores
### Step 6: Heatmaps & Intersections
