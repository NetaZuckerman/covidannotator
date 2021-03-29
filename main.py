import sys, os
import pandas as pd
from datetime import datetime


def main(argv):
    now = datetime.now()
    date_time = now.strftime("%Y%m%d")
    f = open("AfterJoinedList.txt", "r")
    lines = [line.rstrip() for line in f]
    f.close()

    f = open("AfterJoinedList.txt", "a")
    # Search for all the insertions csv files in a given path
    insertionFiles = [os.path.join(root, name)
                      for root, dirs, files in os.walk(argv[0])
                      for name in files
                      if dirs not in lines
                      and "fasta.insertions" in name]
    # Editing the csv files and insert them to an empty list
    li = []
    for table in insertionFiles:
        # Add Files that combined to a list
        for x in table.split('\\'):
            if "NGS" in x or "ngs" in x:
                f.write(x + '\n')
        df = pd.read_csv(table)
        for title in df.columns:
            if "insertion" in title:
                # Remove all the unnecessary title, write only the position of the insertion
                justPosition = title.split()[-1]
                df.rename(columns={title: justPosition}, inplace=True)
                # Change strain to 0 for sorting
            elif "strain" in title:
                df.rename(columns={title: '0'}, inplace=True)
        li.append(df)

    # combine all files in the list
    combined_csv = pd.concat(li, ignore_index=True)
    # sorting columns by their insertion position
    combined_csv = combined_csv[sorted(combined_csv.columns, key=lambda x: tuple(map(int, x.split('_'))))]
    # changing back the strain col name
    combined_csv.rename(columns={'0': 'strain'}, inplace=True)
    # deleting cols that all of their rows are nan
    combined_csv.dropna(axis=1, how='all', inplace=True)
    # save to csv
    combined_csv.to_csv("CombinedCsv_" + date_time + ".csv", index=False, encoding='utf-8-sig')


if __name__ == '__main__':
    main(sys.argv[1:])
