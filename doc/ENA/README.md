# ENA submissions

## Overview

There have been two submissions: 169 and 182. These have been combined under the same umbrella project in submission 183.

## Update

### 2024-08-07

The samples needed to be updated to add now mandatory sample attributes.

The `src/R/UPSC-0182.Sample.R` script was used to regenerate the sample information. A new submission file was created (UPSC-0182.Update.Submission.xml). It was validated and submitted as follows (pwd obviously removed)

```bash
cd doc/ENA
curl -uWebin-763:PWD -k -F "SUBMISSION=@UPSC-0182.Update.Submission.xml" -F "SAMPLE=@UPSC-0182.Updated.Sample.xml" "https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit" > UPSC-0182.Update.Test.Receipt.xml
curl -uWebin-763:PWD -k -F "SUBMISSION=@UPSC-0182.Update.Submission.xml" -F "SAMPLE=@UPSC-0182.Updated.Sample.xml" "https://www.ebi.ac.uk/ena/submit/drop-box/submit" > UPSC-0182.Update.Receipt.xml
```