#### An example of data coming from a memory irradiated with 14-MeV neutrons.
#### It was exposed three times with different patterns. It was read and corrected
#### multiple times, so there is a timestamp that can be used to discriminate MCUs.
###
### Adapt your data following the examples of this template.

##### The memory, CY62167EV, with 16 Mbit, was configured as 2Mx8 instead of
##### 1Mx16. 2 MegaWords require 21 bits. Therefore, the following two constants
##### are easy to understand.
const NbitsAddress=21;
const NWordWidth=8;

##### The CY62167EV is internally divided in two blocks multiplexed by A20, which
##### is the most significant address bit. The following constant indicates that
##### the Nbits4blocks first bits of the address word multiplex separated blocks.
##### If you do not know what to write, define the constant as 0.
const Nbits4blocks = 1;

##### A matrix indicating the written pattern. According to the following matrix,
##### experiment 1 used 0x00 pattern, experiment 2, 0x55, and experiment 3, 0xFF.
##### Respect this notation even if you have performed only one experiment.
const Pattern = [
1 0x00
2 0x55
3 0xFF
];

###############################################################################
### It is possible that you know some critical DV values from previous experiments
### The following two vectors are used to improve the classification of events.
### Respect the way they are noted down. There are some practical examples that
### had been discovered by our research group. Uncomment it.

#### XOR operation
const PreviouslyKnownXOR=[
#0x000002
#0x000006
#0x00c000
#0x00e000
];
#### Positive Subtraction
const PreviouslyKnownPOS=[
#2
#16384
];

###############################################################################

#### The content. The most important part of the date. We tested this memory
#### three times, reading around every minute and correcting errors once discovered.
#### Tests fall into the cathegory of PSEUDOSTATIC.
#### Data coming from experiments require 3 vertical columns. Therefore, columns
#### 1-3 are related to the first experiment (Remember pattern definition above),
#### columns 4-6 to the second one, and 7-9 to the third.
#### Column 1 contains the address where a bitflip was observed.
#### Column 2 contains the read data.
#### Column 3 contains the round when it was observed.
#### That means, for example, that the address 0x1FB1E0 contained an erroneous
#### word, 0x04 instead of pattern=0x00, during Round 1 of experiment 1. Later,
#### 0x002B2D, 0x06E3CC, 0x0A8BFE, and 0x14B49C were discovered in a second round,
#### etc.
####
#### Some important details:
#### 1) as usually it is impossible to obtain exactly the same number of corrupted
####    addresses in the experiments, the addresses columns (1,4, 7, ...) must
####    be completed with the absurd 0xFFFFFFFF to format the matrix. In this
####    example, see columns 4 & 7. The exception is Column 1, in which the number
####    of corrupted addresses was the highest.
#### 2) You may write whatever you want in the DATA and READ-WRITE Cycle columns
####    if a fake value was written in the corresponding address column. Just
####    respect format. See the end of Columns 4-9 in this example.
#### 3) It does not matter how the rounds are sorted. Even if you mix rows (with
####    the notable exception of those containing fake address, which must be
####    always at the column end), the program will work (At least, this is what
####    I want to believe).
#### 4) If you have read only once, you must fill cells of the RWCyC column
####    with 1.
####
Content=round(UInt32, readcsv("DATA/CY62167_090nm_TS.csv"))
