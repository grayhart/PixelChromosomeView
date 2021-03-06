"""
    Copyright 2021 by Neil Millikin and Graham Hart
    contact: neil.millikin@gmail.com graham@the-harts.co.uk
    license: GPLv3

    This program graphically displays DNA match data for three or more siblings

    The comparison of three or more siblings' DNA data has been shown to allow
    derivation of the maternal and paternal grandparents for each sibling
    chromosome pair.

    The process for analyzing similar data sets has been termed "Visual Phasing"

    To use this program,  original raw DNA must be downloaded from
    an original DNA testing vendor, such as AncestryDNA, MyHeritage, 23andMe
    or similar sources.  Once inserted into this program,
    the data is parsed, analyzed and graphically represented.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

""" 
    This version of the application has the following differences 
    from the early version of Neil's application:
        Has only been written for Ancestry DNA raw data files
        Has an option to remove Inserts Deletes and No Calls
        Uses Pandas Dataframes rather than Dictionaries for most of the work
        Creates a graphic line with the colour segments of the matches
            Red - No match
            Yellow - Half match (HIR)
            Green - Full match (FIR)
        Has options to allow configuration of limits for choosing red and green
            How many red or green allele matches constitute a valid segment
            What distance is allowed before a segment is considered to be ended
        Everything is set to Yellow and then Reds are applied and then Greens  
        
"""
import csv
import os
import sys
import inspect
import pprint
import pandas as pd
import numpy as np

from itertools import islice
from collections import OrderedDict
from itertools import combinations

from PIL import Image, ImageDraw, ImageFont

from pixel_config import *

pp = pprint.PrettyPrinter(indent=4)


def lp_2(line_no, name, value):
    if VERBOSITY > 1:
        print("{0}_{1} = {2}\n".format(line_no, name, value))


def lp_3(line_no, name, value):
    if VERBOSITY > 2:
        print("{0}_{1} = {2}\n".format(line_no, name, value))


def pp_2(line_no, description, data_structure):
    if VERBOSITY > 1:
        print("\n{0} -- {1} ==>".format(line_no, description))
        pp.pprint(data_structure)
        print("")


def pp_3(line_no, description, data_structure):
    if VERBOSITY > 2:
        print("\n{0} -- {1} ==>".format(line_no, description))
        pp.pprint(data_structure)
        print("")


# Read in all the raw data files
# Filter by the chromosome we want
# Remove unwanted values
# Only keep common values across the kits
def get_chromsome_dataframe_for_matches(
        all_matches_list):
    this_dir = os.path.dirname(os.path.realpath('__file__'))
    data_dir_name = DATA_FILE_DIRECTORY
    data_file_dir = os.path.join(this_dir, "{0}".format(data_dir_name))
    source_data_file_names = os.listdir(data_file_dir)

    raw_file_names = [f for f in source_data_file_names if 'raw' in f]

    chr_df = pd.DataFrame()
    name_arry = []
    for known_relative in all_matches_list:
        try:
            kr_raw_file_name = [rfn for rfn in raw_file_names
                                if known_relative in rfn
                                and rfn.split('.')[-1] == 'txt'][0]

            lp_2(str(inspect.stack()[0][2]), "Reading kr_raw_file_name raw data",
                 str(kr_raw_file_name))

        except IndexError as e:
            lp_2(str(inspect.stack()[0][2]),
                 "no such file found for {0}".format(known_relative), str(e))
            sys.exit()

        this_kr_raw_file = os.path.join(data_file_dir, kr_raw_file_name)

        # 
        this_df = pd.read_csv(this_kr_raw_file,
                              sep='\t',
                              skip_blank_lines=True,
                              comment='#',
                              header=0,
                              names=['rsid', 'chromosome', 'position',
                                     "{0}-allele1".format(known_relative),
                                     "{0}-allele2".format(known_relative)])

        # Push the columns onto an array so we can remove unique lines later
        name_arry.append("{0}-allele1".format(known_relative))
        name_arry.append("{0}-allele2".format(known_relative))        

        #  Filter out anything except the chromosome we are interested in
        this_df = this_df[this_df['chromosome'] == CHROMOSOME_TO_RENDER]

        # Remove any Inserts or Deletes
        # TODO - This could be changed to only include ACTG which would also remove no calls and 0 
        if FILTER_INSERTS:
            this_df = this_df.drop(this_df[this_df["{0}-allele1".format(known_relative)] == "I"].index)
            this_df = this_df.drop(this_df[this_df["{0}-allele2".format(known_relative)] == 'I'].index)
        if FILTER_DELETES:
            this_df = this_df.drop(this_df[this_df["{0}-allele1".format(known_relative)] == 'D'].index)
            this_df = this_df.drop(this_df[this_df["{0}-allele2".format(known_relative)] == 'D'].index)
        if FILTER_NO_CALLS:
            this_df = this_df.drop(this_df[this_df["{0}-allele1".format(known_relative)] == '0'].index)
            this_df = this_df.drop(this_df[this_df["{0}-allele2".format(known_relative)] == '0'].index)


        # Merge this persons DNA dataframe into a common dataframe keying on Chromosome and Position
        # Only matching chromosome/position pairs will be included in the match
        # NB: Don't use rsid here because even within ancestry DNA kits, the same rposition may be allocated a different rsid
        if len(chr_df.index) > 1:
            if FILTER_MISSING_SNPS:
                chr_df = pd.merge(chr_df, this_df, how='inner', on=('chromosome', 'position'))
            else:
                chr_df = pd.merge(chr_df, this_df, how='outer', on=('chromosome', 'position'))

                # Change NaNs to a character
            chr_df.fillna(NAN_REPLACEMENT_STRING, inplace=True)

        else:
            chr_df = this_df


    if FILTER_COMPLETELY_MATCHED_SEGMENTS:
        print("Removing single allele rows")
        chr_df['allelecount'] = chr_df[name_arry].stack().groupby(level=0).apply(lambda x: x.unique().size)
        chr_df = chr_df[chr_df['allelecount'] > 1]

    lp_2(str(inspect.stack()[0][2]), "Final DF Length: ",
         str(len(chr_df)))

    return chr_df


CHROMOSOME_TO_RENDER = None


def get_match_pair_combinations(
        siblings_to_render,
        extra_match):
    siblings_to_render = list(dict.fromkeys(siblings_to_render))
    all_matches_list = siblings_to_render
    
    match_pair_combinations = [
        (siblings_to_render[m[0]], siblings_to_render[m[1]]) for m in
        list(combinations(range(len(siblings_to_render)), 2))]
    if extra_match and extra_match not in siblings_to_render:
        for sib in siblings_to_render:
            match_pair_combinations.append((sib, extra_match))

        all_matches_list.append(extra_match)

    pp_2(str(inspect.stack()[0][2]), "match_pair_combinations",
         match_pair_combinations)

    return match_pair_combinations, all_matches_list


def insert_combo_match_type_into_common_key_SNP_dict(
        chr_df,
        match_pair_combinations):
    # Loop through each match pair and work out the match type and color
    # and store it in the df

    for match_pair_combination in match_pair_combinations:
        lp_2(str(inspect.stack()[0][2]), "Calculating type and colour of match pair",
             "{0} - {1}".format(match_pair_combination[0], match_pair_combination[1]))

        match1Allele1 = "{0}-allele1".format(match_pair_combination[0])
        match1Allele2 = "{0}-allele2".format(match_pair_combination[0])
        match2Allele1 = "{0}-allele1".format(match_pair_combination[1])
        match2Allele2 = "{0}-allele2".format(match_pair_combination[1])

        # Add a column for this pair 
        newcol = "{0}-{1}".format(match_pair_combination[0],
                                  match_pair_combination[1])

        newColorCol = "{0}-{1}-color".format(match_pair_combination[0],
                                             match_pair_combination[1])
        # When multiple conditions match, the first one is used. So include more specific cases first
        conditions = [
            # Check for any no calls
            (chr_df[match1Allele1] == '0') |
            (chr_df[match1Allele2] == '0') |
            (chr_df[match2Allele1] == '0') |
            (chr_df[match2Allele2] == '0'),

            # Check for any deletes
            (chr_df[match1Allele1] == 'D') |
            (chr_df[match1Allele2] == 'D') |
            (chr_df[match2Allele1] == 'D') |
            (chr_df[match2Allele2] == 'D'),

            # Check for any inserts
            (chr_df[match1Allele1] == 'I') |
            (chr_df[match1Allele2] == 'I') |
            (chr_df[match2Allele1] == 'I') |
            (chr_df[match2Allele2] == 'I'),

            # check for any NaNs (missing SNPs in at least one kit)
            (chr_df[match1Allele1] == NAN_REPLACEMENT_STRING) |
            (chr_df[match1Allele2] == NAN_REPLACEMENT_STRING) |
            (chr_df[match2Allele1] == NAN_REPLACEMENT_STRING) |
            (chr_df[match2Allele2] == NAN_REPLACEMENT_STRING),

            # Chek for fully matched SNPs
            ((chr_df[match1Allele1] == chr_df[match2Allele1]) &
             (chr_df[match1Allele2] == chr_df[match2Allele2])) |
            ((chr_df[match1Allele1] == chr_df[match2Allele2]) &
             (chr_df[match1Allele2] == chr_df[match2Allele1])),

            # Check for no match SNPs
            (chr_df[match1Allele1] != chr_df[match2Allele1]) &
            (chr_df[match1Allele2] != chr_df[match2Allele2]) &
            (chr_df[match1Allele1] != chr_df[match2Allele2]) &
            (chr_df[match1Allele2] != chr_df[match2Allele1])

        ]
        # Add the correct match types to the df
        matchTypeChoices = ['noCallSNP', 'deleteSNP', 'insertSNP', 'missingSNP', 'fullMatchSNP', 'noMatchSNP']
        chr_df[newcol] = np.select(conditions, matchTypeChoices, default='halfMatchSNP')

        # Add the correct colours to the df
        colorChoices = [NO_CALL_SNP_COLOR, DELETE_SNP_COLOR, INSERT_SNP_COLOR, MISSING_SNP_COLOR,
                        FULLY_IDENTICAL_SNP_COLOR, NO_MATCH_SNP_COLOR]
        chr_df[newColorCol] = np.select(conditions, colorChoices, default=HALF_IDENTICAL_SNP_COLOR)

    return chr_df


def calculate_segment_matches(
        chr_df,
        match_pair_combinations):
    # For each matching pair, calculate the strips of green and red (Everything else will be 
    # yellow or gray)
    for match_pair_combination in match_pair_combinations:
        mp_abbr = "{0}-{1}".format(match_pair_combination[0],
                                   match_pair_combination[1])

        mp_abbr_color = "{0}-{1}-color".format(match_pair_combination[0],
                                               match_pair_combination[1])

        mp_abbr_seg_color = "{0}-{1}-seg_color".format(match_pair_combination[0],
                                                       match_pair_combination[1])

        lp_2(str(inspect.stack()[0][2]), "Creating segments for ", str(mp_abbr))

        # Initialise a new column with yellows as default
        chr_df[mp_abbr_seg_color] = 'yellow'

        # Create a df with just the fields we need so it isn't so large
        red_df = chr_df[chr_df[mp_abbr_color] == 'red']
        green_df = chr_df[chr_df[mp_abbr_color] == 'limegreen']

        # Loop through the reds record a start point and stop when the gap to the next one is 
        # bigger than the RED_GAP
        start = -1
        current_position = -1
        count = 0
        for index, row in red_df.iterrows():
            if start == -1:
                # This entry is the start of a new segment
                if index < RED_GAP:
                    start = index
                    count += 1
                else:
                    start = index
                    count += 1
                current_position = index
            else:
                # set up some values 
                if index - current_position > RED_GAP:
                    # We are at the end of a segment of reds set up to the previous setting
                    if count >= MIN_REDS:
                        chr_df = set_segmet_red(start, current_position,
                                                mp_abbr_seg_color, 'red', chr_df)
                    start = row.name
                    current_position = index
                    count = 0
                else:
                    count += 1
                    current_position = index
        # If we're at the end then see if we have a last segment to paint
        if count >= MIN_REDS:
            chr_df = set_segmet_red(start, current_position,
                                    mp_abbr_seg_color, 'red', chr_df)

        # Loop through the greens record a start point and stop when the gap to the next one is
        # bigger than the GREEN_GAP
        start = -1
        current_position = -1
        for index, row in green_df.iterrows():
            if start == -1:
                # This entry is the start of a new segment
                if index < GREEN_GAP:
                    start = index
                    count += 1
                else:
                    start = index
                    count += 1
                current_position = index
            else:
                # set up some values 
                if index - current_position > GREEN_GAP:
                    # We are at the end of a segment of reds set up to the previous setting
                    if count >= MIN_GREENS:
                        chr_df = set_segmet_red(start, current_position,
                                                mp_abbr_seg_color, 'limegreen', chr_df)
                    start = row.name
                    current_position = index
                    count = 0
                else:
                    count += 1
                    current_position = index

        # If we're at the end then see if we have a last segment to paint
        if count >= MIN_GREENS:
            chr_df = set_segmet_red(start, current_position,
                                    mp_abbr_seg_color, 'limegreen', chr_df)

        # We need to decide if we colour the last one
    return chr_df


def set_segmet_red(start_position,
                   end_position,
                   mp_abbr_seg_color,
                   color,
                   chr_df):
    # Set all the segments between start and end position to red for this match pair
#    chr_df.loc[start_position:end_position][mp_abbr_seg_color] = color
    chr_df.loc[start_position:end_position, mp_abbr_seg_color] = color
    
#    dfmi.loc[:, ('one', 'second')]
    
    return chr_df


def create_comparison_base_strip_image(
        width,
        height):
    comarison_base_strip_image = Image.new(
        'RGB',
        (width, height),
        color='white')
    return comarison_base_strip_image


def create_chrom_whole_page_image(
        width,
        height,
        color):
    chromosome_full_page_image = Image.new(
        'RGB',
        (width, height),
        color=color)
    return chromosome_full_page_image


def draw_single_SNP_line(
        width,
        height,
        color):
    single_SNP_line = Image.new(
        'RGB',
        (width, height),
        color=color)
    return single_SNP_line


def show_match_graphics(
        chr_df,
        match_pair_combinations):
    file_lines = len(chr_df.index)

    title = "Pixel View GEDmatch Chr {0}"
    file_prefix = ""
    if FILTER_INSERTS:
        title += " - no Inserts"
        file_prefix += "no I"
    else:
        title += " - Inserts"
        file_prefix += "-I"

    if FILTER_DELETES:
        title += " - no Deletes"
        file_prefix += "-no D"
    else:
        title += " - Deletes"
        file_prefix += "-D"

    if FILTER_NO_CALLS:
        title += " - no NoCalls"
        file_prefix += "-no NoCalls"
    else:
        title += " - NoCalls"
        file_prefix += "-NoCalls"

    if FILTER_COMPLETELY_MATCHED_SEGMENTS:
        title += " - no CompletelyMatchedSNPs"
        file_prefix += "-No CompleteMatched"
    else:
        title += " - CompletelyMatchedSNPs"
        file_prefix += "-CompletelyMatched"

    matches_to_show = len(match_pair_combinations)
    chrom_whole_page_image = create_chrom_whole_page_image(

        max((file_lines
             + CHROM_PAGE_LEFT_BORDER
             + CHROM_PAGE_RIGHT_BORDER),
            MINIMUM_PAGE_WIDTH),

        (CHROM_PAGE_TOP_BORDER
         + CHROM_PAGE_BOTTOM_BORDER
         + CHROM_PAGE_TITLE_SPACE)
        + (matches_to_show * 2 * (
                HEIGHT_OF_CHROMOSOME_IMAGE
                + SPACE_BETWEEN_MATCHES)),

        BACKGROUND_COLOR)

    lp_3(str(inspect.stack()[0][2]), "matches_to_show", str(matches_to_show))

    page_draw = ImageDraw.Draw(chrom_whole_page_image)
    page_draw.text((CHROM_PAGE_LEFT_BORDER, CHROM_PAGE_TOP_BORDER), title,
                   font=ImageFont.truetype("arial.ttf",
                                           CHROM_PAGE_TEXT_FONT_SIZE), fill='black')

    match_shown_number = 0
    
    
    # Create a dataframe with only the fields we need
    col_list = []
    for match_pair_combination in match_pair_combinations:
        mp_abbr = "{0}-{1}".format(match_pair_combination[0],
                                   match_pair_combination[1])

        mp_abbr_color = "{0}-{1}-color".format(match_pair_combination[0],
                                               match_pair_combination[1])

        mp_abbr_seg_color = "{0}-{1}-seg_color".format(match_pair_combination[0],
                                                       match_pair_combination[1])
        col_list.append(mp_abbr)
        col_list.append(mp_abbr_color)
        col_list.append(mp_abbr_seg_color)
        

    # only keep the col
    chr_df = chr_df[col_list]

    for match_pair_combination in match_pair_combinations:
        mp_abbr = "{0}-{1}".format(match_pair_combination[0],
                                   match_pair_combination[1])

        mp_abbr_color = "{0}-{1}-color".format(match_pair_combination[0],
                                               match_pair_combination[1])

        mp_abbr_seg_color = "{0}-{1}-seg_color".format(match_pair_combination[0],
                                                       match_pair_combination[1])

        lp_2(str(inspect.stack()[0][2]), "Creating graphics for ", str(mp_abbr))

        # Create a df with just the fields we need so it isn't so large
        this_df = chr_df[[mp_abbr_color, mp_abbr_seg_color]]

        base_position = 0

        comparison_base_strip_image = create_comparison_base_strip_image(
            file_lines, HEIGHT_OF_CHROMOSOME_IMAGE)

        segment_base_strip_image = create_comparison_base_strip_image(
            file_lines, HEIGHT_OF_CHROMOSOME_IMAGE)

        def chr_df_iter(comparison_base_strip_image,
                        match_color,
                        seg_color,
                        base_position
                        ):
            color = match_color
            seg_color = seg_color

            paste_position = (base_position, 0)

            single_pixel_line = draw_single_SNP_line(WIDTH_OF_SNP_LINE,
                                                     HEIGHT_OF_CHROMOSOME_IMAGE, color)

            single_pixel_seg_line = draw_single_SNP_line(WIDTH_OF_SNP_LINE,
                                                         HEIGHT_OF_CHROMOSOME_IMAGE, seg_color)

            comparison_base_strip_image.paste(single_pixel_line,
                                              (paste_position))

            segment_base_strip_image.paste(single_pixel_seg_line,
                                           (paste_position))

            chr_df.at[base_position, 'color'] = color
            base_position += 1

        this_df.apply(lambda row, base_position=base_position: chr_df_iter(comparison_base_strip_image,
                                                                           row[mp_abbr_color],
                                                                           row[mp_abbr_seg_color],
                                                                           row.name), axis=1)

        chrom_whole_page_image.paste(
            comparison_base_strip_image,

            (CHROM_PAGE_LEFT_BORDER,

             (CHROM_PAGE_TOP_BORDER + CHROM_PAGE_TITLE_SPACE
              + (2 * match_shown_number * (
                             HEIGHT_OF_CHROMOSOME_IMAGE
                             + SPACE_BETWEEN_MATCHES)))))

        chrom_whole_page_image.paste(
            segment_base_strip_image,

            (CHROM_PAGE_LEFT_BORDER,

             (CHROM_PAGE_TOP_BORDER + CHROM_PAGE_TITLE_SPACE
              + ((2 * match_shown_number + 1) * (
                             HEIGHT_OF_CHROMOSOME_IMAGE
                             + SPACE_BETWEEN_MATCHES))
              - SPACE_BETWEEN_MATCHES)))

        draw = ImageDraw.Draw(chrom_whole_page_image)

        arial = ImageFont.truetype("arial.ttf", CHROM_PAGE_TEXT_FONT_SIZE)

        draw.text((
            CHROM_PAGE_TEXT_BORDER,

            (CHROM_PAGE_TOP_BORDER
              + CHROM_PAGE_TITLE_SPACE
                + (2 * match_shown_number
                  * (HEIGHT_OF_CHROMOSOME_IMAGE
                   + SPACE_BETWEEN_MATCHES)))),            
            mp_abbr, font=arial, fill='black')

        match_shown_number += 1

    this_dir = os.path.dirname(os.path.realpath('__file__'))
    image_file_dir = os.path.join(this_dir, "{0}".format(IMAGE_FILE_DIRECTORY))

    image_file = "{0}\\pixel_chr_{3}_{2}_chr_{1}.png".format(image_file_dir,
                                                  CHROMOSOME_TO_RENDER,
                                                  file_prefix,
                                                  '-'.join(siblings_to_render))
    print("Saving image")
    chrom_whole_page_image.save(image_file)

    lp_2(str(inspect.stack()[0][2]), "Finished graphics for chromosome: ", str(CHROMOSOME_TO_RENDER))
    print('Image saved as:{0}'.format(image_file))


if __name__ == '__main__':

    def get_valid_chromosome_number():
        global CHROMOSOME_TO_RENDER
        SUGGESTED_CHROMOSOME_TO_RENDER = input("""
            "Enter a valid chromosome number to process 1-22  
            
            or type 'quit' to exit program     """)

        if SUGGESTED_CHROMOSOME_TO_RENDER == 'quit':
            print("goodbye")
            sys.exit()

        elif SUGGESTED_CHROMOSOME_TO_RENDER in [str(i) for i in range(1, 23)]:
            CHROMOSOME_TO_RENDER = int(SUGGESTED_CHROMOSOME_TO_RENDER)

        else:
            get_valid_chromosome_number()


    if SINGLE_CHROMOSOME:
        get_valid_chromosome_number()
        start_chr = CHROMOSOME_TO_RENDER
        end_chr = CHROMOSOME_TO_RENDER + 1
        lp_2(str(inspect.stack()[0][2]), "Option for single chromosome selected. Chromosome: ", str(CHROMOSOME_TO_RENDER))

    else:
        start_chr = 1
        end_chr = 24
        lp_2(str(inspect.stack()[0][2]), "Option for all chromosomes selected.", '')

    match_pair_combinations, all_matches_list = get_match_pair_combinations(siblings_to_render, extra_match)
 
    for CHROMOSOME_TO_RENDER in range(start_chr, end_chr):
        lp_2(str(inspect.stack()[0][2]), "Processing chromosome: ", str(CHROMOSOME_TO_RENDER))


        chr_df = None
        chr_df = get_chromsome_dataframe_for_matches(all_matches_list)

        # Reset the index so it goes from 0 on
        chr_df.reset_index(inplace=True, drop=True)

        chr_df = insert_combo_match_type_into_common_key_SNP_dict(
            chr_df,
            match_pair_combinations)

        chr_df = calculate_segment_matches(
            chr_df,
            match_pair_combinations)

        show_match_graphics(
            chr_df,
            match_pair_combinations)
