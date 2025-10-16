# Deduper Strategy & Pseudocode

### **Defining the Problem**
Given a `.sam` (single-end) file of uniquely mapped reads, and a `.txt` file containing the known UMIs (unique molecular identifiers), we want to strategically remove all PCR duplicates and retain ONLY the first read encountered if duplicates are found. Additionally, we want to **avoid** overloading our memory with millions of reads while processing. Before approaching this problem, let's understand what reads are considered duplicates.

#### Identifying Duplicates
We will utilize the `.sam` file to help us identify duplicates. Duplicates are identified if:
- The alignment position is the same:
	- Same chromosome number (`sam file column 3`)
	- Same strand (`sam file column 2: Bitwise Flag 4 [Mapped or Unmapped] & Flag 16 [+ or -]`)
	- Same 5' starting position (1-based leftmost mapping POSition (`sam file column 4`) + CIGAR string (`sam column 6`) if there's soft clipping) (look into soft clipping below)
		- Soft clipping:
			- The CIGAR string at `sam file column 6` will tell us if there's soft clipping with the notation `S`. This, alongside the strandness information, will help us determine the true 5' starting position!
- The Same UMI!

### Key Considerations/ Strategy
- Sort the sam file by their chromosome number first. Each chromosome number will be tracked with a `profile` that is a `dict` with format `{"UMI": str, "strand": +/-,"5_prime": int}`.  This will be put inside a `set`. 
	- This `set` will  track for duplicates as we parse each record in that same chromosome number. If we encounter a unique record, the list will append a new profile and will write that record to the output file.
	- Once we're done parsing through a chromosome and ready for the next one, the `set` of `profile` will be cleared out and reset to prevent memory overload.
- If an unknown UMI is encountered, disregard that record and proceed to the next
- If the read is unmapped, disregard that record. <- Not Sure if we need this
- Retain only the first copy if a duplicate is encountered

### The Pseudocode

1. Sort the input sam file by their chromosome number.
2. Argument Parsing
	1. The function `get_args()` will take in the following `argparse` options:
		1. `-f`: absolute path to sorted sam file
		2. `-o`: absolute path to deduplicated sam file
		3. `-u`: path to the file containing the list of UMIs
		4. `-h`: prints useful help message of the script
	2. Store the passed arguments into the following variables:
		1. `input_file <- args.f`
		2. `output_file <- args.o`
		3. `umi_file <- args.u`
3. Call the function `parse_umi_txt()` to generate a list of the valid UMIs and save that to the variable `umi_list`
4. Create a global variable called `profiles` which is a `set{}` of dictionaries of format `{"UMI": str, "strand": +/-,"5_prime": int}`.
5. Create a global variable called `chr_track` and set it to `None` to keep track of the chromosome number we're dealing with.
6. Open the `input_file` and `output_file` and use a `for-loop` to iterate line by line, parsing each record. Check for the following step:
	1. Step 1a: Validate the umi by calling the function `determine_umi()`. If null, skip to the next record. If valid, store it to variable `umi` and proceed to the next step.
	2. Step 1b: Validate if the record is mapped or not by calling `mapped()`. If mapped, proceed to the next step. Else, skip to the next record.
	3. Step 2: Retrieve the chromosome number by calling `chr_num`. Store it to the variable `chr`.
		1. If `chr != chr_track`, set `chr_track = chr` and clear `profiles`
	4. Step 3: Retrieve the strandness by calling `direction()` and store it to variable `strand`.
	5. Step 4: Determine the 5' starting position by calling `five_prime_pos` and store it to the variable `true_start`
	6. Step 5: Store the variables `umi, strand, true_start` as a dict  called `info` `{"UMI": umi, "strand": strand,"5_prime": true_start}`
	7. Step 6: Check for duplication and write records:
		1. If `profiles` is empty, append `info` to it and write the record to the output file
		2. Else, check if `info` is in the set. If yes, skip to the next record.
	8. Step 7: Repeat Step 1-6 until EOF.
7. Ta-da!
### Defining High Level Functions
Functions that involve parsing a `sam_file` record will use `string` functions `split`. For functions that involve looking into the CIGAR string, we will use `regex` to help split at every alphabet.
#### Parse UMI 
```python
def parse_umi_txt(umi_file: str):
	"""
	Take the UMI txt file and return a list of the 96 UMIs in strings
	"""
Input: Umi.txt
Output: ["ACTGTT", ..."CCTGTA"]
```
#### Determine UMI
```python
def determine_umi(sam_record: str, umi_list: list[str]):
	"""
	Take a sam record file and the UMI list and return the UMI number only if it's valid. Returns None if invalid.
	"""
Input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC... ["CTGTTCAC"...]
Output: CTGTTCAC

Input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAA... ["CTGTTCAC"...]
Output: None
```
#### Determine Mapped Reads
```python
def mapping(sam_record: str):
	"""
	Takes a single sam_record and return True if mapped, else false.
	"""
Input: QNAME 104 7 123456 36 71M * 0 0 ACTGTS...TGS ...
Expected: False

Input: QNAME 100 7 123456 36 71M * 0 0 ACTGTS...TGS ...
Expected: True
```
#### Determine Chromosome Number
```python
def chr_num(sam_record: str):
	"""
	Takes a single sam_record and return the chromosome number (3rd col).
	"""
Input: QNAME 100 7 123456 36 71M * 0 0 ACTGTS...TGS ...
Expected: 7
```
#### Determine Left-Most Mapping Position
```python
def left_pos(sam_record: str):
	"""
	Takes a single sam_record and return the 1-based leftmost mapping position (4th col).
	"""
Input: QNAME 100 7 123456 36 71M * 0 0 ACTGTS...TGS ...
Expected: 123456
```
#### Determine Strandness
```python
def direction(sam_record: str):
	"""
	Takes a single sam_record and return the strandness number (3rd col). Will use bitwise operator "&".
	"""
Input: QNAME 100 7 123456 36 71M * 0 0 ACTGTS...TGS ...
Expected: "+"

Input: QNAME 116 7 123456 36 71M * 0 0 ACTGTS...TGS ...
Expected: "-"
```
#### Determine Right-Most Mapping Position
```python
def right_pos(sam_record: str):
	"""
	Takes a single sam_record and return the 1-based rightmost mapping position. Will parse the CIGAR string and call left_pos(). Strand doesn't matter here. Will use regex expression to split the CIGAR to ignore soft-clips notation focusing only on the values in between to determine the length.
	"""
# Example for + strand
Input: QNAME 100 7 123456 36 71M * 0 0 ACTGTS...TGS ...
Expected: 123527

# Example for - strand
Input: QNAME 116 7 123456 36 20M30D21M * 0 0 ACTGTS...TGS ...
Expected: 123527

# Example for + strand with soft clipping 
Input: QNAME 100 7 123456 36 5S45D15M6S * 0 0 ACTGTS...TGS ...
Expected: 123516

# Example for - strand with soft clipping at the end
Input: QNAME 116 7 123456 36 5S45D15M6S * 0 0 ACTGTS...TGS ...
Expected: 123516
```
#### Adjusting for SoftClipping
```python
def adj_soft_clip(sam_record: str):
	"""
	Takes a single sam_record and adjust for soft clipping if needed. Will parse the CIGAR string and return an int value. Will call the function direction() to determine if necessary.
	"""
# Example for + strand with no soft-clipping
Input: QNAME 100 7 123456 36 69M2S * 0 0 ACTGTS...TGS ...
Expected: 0 # default value for no soft clipping

# Example for + strand with soft-clipping
Input: QNAME 100 7 123456 36 10S71M12S * 0 0 ACTGTS...TGS ...
Expected: 10

# Example for - strand with no soft-clipping at the end
Input: QNAME 116 7 123456 36 2S69M * 0 0 ACTGTS...TGS ...
Expected: 0

# Example for - strand with soft-clipping at the end
Input: QNAME 116 7 123456 36 10S71M12S * 0 0 ACTGTS...TGS ...
Expected: 12
```
#### Determine 5' start pos
```python
def five_prime_pos(sam_record: str):
	"""
	Takes a single sam_record and return the 5' starting position. This function calls the adj_soft_clip() and right_pos(). If no soft clip, just call right_pos!
	"""
# Example for + strand with no soft-clipping
Input: QNAME 100 7 123456 36 71M * 0 0 ACTGTS...TGS ...
Expected: 123456

# Example for + strand with soft-clipping
Input: QNAME 100 7 123456 36 2S71M * 0 0 ACTGTS...TGS ...
Expected: 123454

# Example for - strand with no soft clipping at the end
Input: QNAME 116 7 123456 36 5S45D21M * 0 0 ACTGTS...TGS ...
Expected: 123522

# Example for - strand with soft clipping at the end
Input: QNAME 116 7 123456 36 5S45D15M6S * 0 0 ACTGTS...TGS ...
Expected: 123522
```
