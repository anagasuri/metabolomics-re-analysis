import pandas as pd

POS_FILE = "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/3. Frequency/POS/CLEANED_final_annotated_POS.csv"
NEG_FILE = "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/3. Frequency/NEG/CLEANED_final_annotated_neg.csv"
OUT_FILE = "/Users/amritanagasuri/Desktop/new Q-TOF Analysis/3. Frequency/CLEANED_final_annotated_exogenous.csv"

ENDO_COL = "If endogenous metabolite"

def has_N_token(series):
    # Uppercase, strip spaces, normalize delimiters (treat comma like semicolon)
    s = (series.astype(str)
         .str.upper()
         .str.replace(r"\s+", "", regex=True)
         .str.replace(",", ";", regex=False))
    # Match 'N' as a token at start, end, or between delimiters ; or ,
    # (after normalization, only semicolons remain)
    return s.str.contains(r'(?:^|;)N(?:$|;)', regex=True, na=False)

def main():
    pos = pd.read_csv(POS_FILE, low_memory=False)
    neg = pd.read_csv(NEG_FILE, low_memory=False)

    pos_N = pos[has_N_token(pos[ENDO_COL])].copy()
    neg_N = neg[has_N_token(neg[ENDO_COL])].copy()

    stacked = pd.concat([pos_N, neg_N], ignore_index=True)
    merged = stacked.drop_duplicates(subset=["Alignment ID"], keep="first")
    merged.to_csv(OUT_FILE, index=False)

if __name__ == "__main__":
    main()
