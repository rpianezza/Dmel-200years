#!/bin/bash

# Define full paths for clarity
popte2="/Volumes/Element/Backup/dmel_twocenturies_local/Revision/final_USA_UK/popte2-v1.10.03.jar"
ppileup="/Volumes/Element/Backup/dmel_twocenturies_local/Revision/final_USA_UK/output_final/USA.ppileup.gz"
output="/Volumes/Element/Backup/dmel_twocenturies_local/Revision/final_USA_UK/output_sh_3"
refgenome="/Volumes/Element/Backup/dmel_twocenturies_local/Supp_fig/PopTE2/Dmel_TE_merged.fa"

# Commands
java -jar "$popte2" identifySignatures \
  --ppileup "$ppileup" \
  --mode separate \
  --output "$output/USA-min-count3.signatures" \
  --min-count 3

java -jar "$popte2" frequency \
  --ppileup "$ppileup" \
  --signature "$output/USA-min-count3.signatures" \
  --output "$output/USA-min-count3.freqsig"

sed 's/,//g' "$output/USA-min-count3.freqsig" > "$output/USA-min-count3_fixed.freqsig"

java -jar "$popte2" filterSignatures \
  --input "$output/USA-min-count3_fixed.freqsig" \
  --output "$output/USA-min-count3.filter.freqsig" \
  --max-otherte-count 2 \
  --max-structvar-count 2

sed 's/,//g' "$output/USA-min-count3.filter.freqsig" > "$output/USA-min-count3_fixed.filter.freqsig"

java -jar "$popte2" pairupSignatures \
  --signature "$output/USA-min-count3_fixed.filter.freqsig" \
  --ref-genome "$refgenome" \
  --hier names_TEs.hier \
  --min-distance -200 \
  --max-distance 300 \
  --output "$output/USA-min-count3.teinsertions"

sed 's/\([0-9]\),\([0-9]\)/\1.\2/g' "$output/USA-min-count3.teinsertions" > "$output/USA-min-count3.teinsertions_fixed"