import pandas as pd
import argparse

color_hash = {
    0 :"229,229,229",
    1 :"196,196,196",
    2 :"0,0,0",
    3 :"0,0,205",
    4 :"65,105,225",
    5 :"100,149,237",
    6 :"180,238,180",
    7 :"255,255,0",
    8 :"255,165,0",
    9 :"139,26,26",
    10 :"255,0,0"
}

def get_color(int):
    return color_hash.get(int, "255,0,0")

def main():
    p = argparse.ArgumentParser()
    p.add_argument('bed')
    p.add_argument('bed9')
    args = p.parse_args()

    # Reads bedGraph file
    df = pd.read_csv(args.bed, delimiter = '\t', names = ['chrom','chromStart','chromEnd','label'])

    # Creates new columns
    df['label'] = [round(x) for x in df['label']]
    df['score'] = pd.Series([0]*len(df))
    df['strand'] = pd.Series(['+']*len(df))
    df['thickStart'] = df['chromStart']
    df['thickEnd'] = df['chromEnd']
    df['itemRgb'] = [get_color(x) for x in df['label']]

    # Saves output as bed 9 format
    df.to_csv(args.bed9, sep = '\t', header = False, index = False)

if __name__ == '__main__':
    main()
