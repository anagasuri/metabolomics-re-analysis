# Modified by Amrita

# for pos and neg clean
# for pos and neg corrected data with reps 

# re-run/modifed input files by AMRITA NAGASURI for re-analysis 9/24/25

import pandas as pd
import re

def remove_rows_with_words(input_file, output_file, words_to_remove, column=None):
    """
    Removes rows (and some columns) from a CSV if they contain any of the given words.

    Rules:
    - If a column's NAME contains a word, drop that column entirely.
    - For the remaining columns, if any cell contains a word, drop that row.
    - If `column` is provided, only that column is used for row filtering (columns are still
      dropped based on their header names).
    """
    # Load CSV as strings
    df = pd.read_csv(input_file, dtype=str).fillna("")

    # Build one regex that matches any of the given words literally (substring match)
    pattern = re.compile("|".join(re.escape(w) for w in words_to_remove))

    # 1) Drop columns whose *header* contains any of the words
    cols_to_drop = [c for c in df.columns if pattern.search(c)]
    if cols_to_drop:
        df = df.drop(columns=cols_to_drop)

    # 2) Drop rows:
    #    - If `column` is specified and present, only look there.
    #    - Otherwise, look across ALL remaining columns.
    if column is not None and column in df.columns:
        row_has_word = df[column].str.contains(pattern, na=False)
    else:
        row_has_word = df.apply(lambda s: s.str.contains(pattern, na=False)).any(axis=1)

    cleaned_df = df.loc[~row_has_word]

    # Save
    cleaned_df.to_csv(output_file, index=False)
    print(f"Cleaned file saved as: {output_file}")


# Example usage
if __name__ == "__main__":
    input_csv = "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/3. Frequency/POS/final_annotated_POS.csv"
    output_csv = "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/3. Frequency/POS/CLEANED_final_annotated_POS.csv"
    words = ["P_BI00707", "P_BI00709",  
             "P_BI00712", "P_BI00715", 
             "P_BI00724", "P_BI00727", 
             "P_BI00730", "P_BI00741", 
             "P_BI00759", "P_BI00746", 
             "P_BI00753", "C_BI00707", "C_BI00709",  
             "C_BI00712", "C_BI00715", 
             "C_BI00724", "C_BI00727", 
             "C_BI00730", "C_BI00741", 
             "C_BI00759", "C_BI00746", 
             "C_BI00753"]

    # Use this if you want only the Spectrum column to control row filtering:
    # remove_rows_with_words(input_csv, output_csv, words, column="Spectrum reference file name")

    # Or this to scan all remaining columns for row filtering:
    remove_rows_with_words(input_csv, output_csv, words)