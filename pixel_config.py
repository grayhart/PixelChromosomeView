# There are THREE MANDATORY settings:

"""
There are TWO MANDATORY settings:
"""

DATA_FILE_DIRECTORY = 'raw_files'
IMAGE_FILE_DIRECTORY = 'images'
siblings_to_render = ['KEN', 'ALAN', 'BRIAN']

# optionally you can compare the siblings with one additional relative
extra_match = 'BERYL'

"""
file names must:
1) contain name of source vendor (i.e. Ancestry, 23andMe, MyHeritage, etc
2) contain person's name exactly as listed in 'siblings_to_render' below
3) contain the word 'raw'
4) lastly, raw data musst must be in  .txt files

example: 23andMe_JULIE_raw_dna.txt

copy all relevant raw dna files into the directory your specify
"""

"""
    Set the SINGLE_CHROMOSOME TO True to be prompted to choose 
    the chromosome you want to process.
    False means it will process all chromosomes
"""
SINGLE_CHROMOSOME = True

# Filter out any entries with 'I', 'D' or '0' ()
FILTER_INSERTS = True
FILTER_DELETES = True
FILTER_NO_CALLS = True
FILTER_MISSING_SNPS = True

# removes a lot of 'noise' SNPs that don't contribute to the analysis
FILTER_COMPLETELY_MATCHED_SEGMENTS = True
COLOUR_COMPLETELY_MATCHED_SEGMENTS = False

# optional settings to change appearance of rendered chromosome pairs  #
FULLY_IDENTICAL_SNP_COLOR = 'limegreen'
NO_MATCH_SNP_COLOR = 'red'
HALF_IDENTICAL_SNP_COLOR = 'yellow'

FULL_MATCH_SEGMENT_COLOR = 'limegreen'
HALF_MATCH_SEGMENT_COLOR = 'yellow'
NO_MATCH_SEGMENT_COLOR = 'red'


DELETE_SNP_COLOR = 'gray'
INSERT_SNP_COLOR = 'gray'
NO_CALL_SNP_COLOR = 'gray'
MISSING_SNP_COLOR = 'orange'
COMPLETELY_MATCHED_SNP_COLOR = 'gray'

BACKGROUND_COLOR = 'lightsteelblue'
CHROMOSOME_BASE_COLOR = 'white'

# Variables to control the segment colouring
RED_GAP = 600
GREEN_GAP = 1

MIN_REDS = 5
MIN_GREENS = 100

WIDTH_OF_SNP_LINE = 1
HEIGHT_OF_CHROMOSOME_IMAGE = 600
SPACE_BETWEEN_MATCHES = 300
CHROM_PAGE_LEFT_BORDER = 2500
CHROM_PAGE_TEXT_BORDER = 20
CHROM_PAGE_TEXT_FONT_SIZE = 360
CHROM_PAGE_RIGHT_BORDER = 50
CHROM_PAGE_TOP_BORDER = 150
CHROM_PAGE_TITLE_SPACE = 600
CHROM_PAGE_BOTTOM_BORDER = 150
MINIMUM_PAGE_WIDTH = 1750

VERBOSITY = 2
'''
    DON'T CHANGE OPTIONS BELOW HERE 
'''

# Replace NaNs with a specific value so I can process them more easily
NAN_REPLACEMENT_STRING = 'Z'
