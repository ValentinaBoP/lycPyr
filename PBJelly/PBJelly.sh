#! /bin/bash -l
#SBATCH -A snic2018-8-266
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 20:00:00
#SBATCH -J PBJELLY_lycPyr5.1
#SBATCH -o PBJELLY_lycPyr5.1.output
#SBATCH -e PBJELLY_lycPyr5.1.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load PBSuite/15.8.24

REF=/home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/2.GapFilling/Data/lycPyr5.1.fasta
QUAL=/home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/2.GapFilling/Data/lycPyr5.1.qual
PROTOCOL=/home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/2.GapFilling/Code/Protocol_lycPyr5.1.xml

#set -e

#fakeQuals.py $REF $QUAL

# this is good is we run the default parameters
#for STAGE in setup mapping support extraction assembly output
#do
#    Jelly.py $STAGE $PROTOCOL
#done

#Jelly.py setup -x --minGap=10 $PROTOCOL

#for STAGE in mapping support extraction

for STAGE in support extraction
do
        Jelly.py $STAGE $PROTOCOL
done

Jelly.py assembly -x --maxWiggle=5000 $PROTOCOL

Jelly.py output $PROTOCOL
