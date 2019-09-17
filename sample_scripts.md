## Sample scripts

#### 1. `minimap2.sh` runs reference indexing and aligning target's sequence to human

It required to download human fasta as a reference.

```bash
#!/bin/bash

cd DATA
wget -nc https://biobox-info.github.io/test_data/MT-human.fa
cd ..

human=MT-human.fa

# 1. build minimap2 index of MT-human
singularity exec -B DATA:/scratch Minimap2_latest.sif \
minimap2 -x map-ont -d /scratch/${human%.fa}-ont.mmi /scratch/$human

```

To download the test data:

```bash
#!/bin/bash

cd DATA
wget -nc https://biobox-info.github.io/test_data/MT-orang.fa
wget -nc https://biobox-info.github.io/test_data/MT-mouse.fa
wget -nc https://biobox-info.github.io/test_data/MT-whale.fa
cd ..

```

The below script is an actual script to run on the multiple input data

```bash
#!/bin/bash

# 2. align & map Target to MT-human
singularity exec -B DATA:/scratch Minimap2_latest.sif \
minimap2 -a /scratch/MT-human-ont.mmi /scratch/$1 -o /scratch/${1%.fa}.sam

# 3. Target sam to bam
singularity exec -B DATA:/scratch Samtools_latest.sif \
samtools view -bhS /scratch/${1%.fa}.sam -o /scratch/${1%.fa}.bam

# 4. create an index Target bam file
singularity exec -B DATA:/scratch Samtools_latest.sif \
samtools index /scratch/${1%.fa}.bam

```

To check the output in IGV:

```python
# 1st Cell
import igv
b = igv.Browser(
    {"reference": {
        "id": "mt",
        "fastaURL": "DATA/MT-human.fa",
        "indexed": False,
    }}
)
b.show()

# 2nd Cell
b.load_track(
    {
        "name": "MT-Orang",
        "url": "DATA/MT-orang.bam",
        "indexURL": "DATA/MT-orang.bam.bai",
        "format": "bam",
        "type": "alignment",
        "height": 100
    })

b.load_track(
    {
        "name": "MT-Mouse",
        "url": "DATA/MT-Mouse.bam",
        "indexURL": "DATA/MT-Mouse.bam.bai",
        "format": "bam",
        "type": "alignment",
        "height": 100
    })
b.load_track(
    {
        "name": "MT-Whale",
        "url": "DATA/MT-Whale.bam",
        "indexURL": "DATA/MT-Whale.bam.bai",
        "format": "bam",
        "type": "alignment",
        "height": 100
    })

```

#### 2. `getSampleData.sh` downloads sample data part from Trinity repo.

```bash
#!/bin/bash -ve

mkdir trinityrnaseq
cd trinityrnaseq
git init
git remote add -f origin https://github.com/trinityrnaseq/trinityrnaseq.git
git config core.sparseCheckout true
echo "sample_data/test_Trinity_Assembly" >> .git/info/sparse-checkout
git pull origin master
mv sample_data/test_Trinity_Assembly/* .
rm -rf sample_data
git remote rm origin
cd ..
mv trinityrnaseq sample_data
```
