% (A) (M) Introduction and take in file name and store contents into variable
fprintf("\nHello! Welcome to Wilton's BIEN101 Project 1!\n")
filename = input("Please input the file name below \n(If using my file it should be: ' wiltondna.txt ') \n(Also make sure it is in MATLAB folder)\n File Name: ", "s") ;
file = fopen(filename,"r") ;
dnatext = fileread(filename) ;
fclose(file) ;

% (B) Take in input for Template or Coding strand, heck input and store file contents in corresponding variables
strandType = upper(input("Please input whether strand is Template or Coding (Enter T or C): ", "s")) ;
if strandType == "T"
    strandType = "Template" ;
    templateStrand = dnatext ; 
    codingStrand = reverse(dnatext) ; % reverse because coding is reverse of template
elseif strandType == "C" 
    strandType = "Coding" ; 
    codingStrand = dnatext ; 
    templateStrand = reverse(dnatext) ; 
end 

% (D) Initialize variables to count bases and their complement amounts to be used later on 

adenineCount = 0 ; 
thymineComplement = 0 ; 

thymineCount = 0 ; 
adenineComplement = 0 ; 

guanineCount = 0 ; 
cytosineComplement = 0 ; 

cytosineCount = 0 ; 
guanineComplement = 0 ; 

% (C) Initialize complementary strand as empty char to add complementary base to
complementaryStrand = ''; 

% Create for-loop to parse and iterate through file contents (stored in "dnatext"), check each element (base) to count and add complement to complementary strand 
for i = 1:length(dnatext)
    n = dnatext(i) ;
    if isletter(n)
        n = upper(n) ; 
        if n == 'A'
            adenineCount = adenineCount + 1 ; 
            thymineComplement = thymineComplement + 1 ; 
            complementaryStrand = [complementaryStrand, 'T'] ; 
        elseif n =='T' 
            thymineCount = thymineCount + 1 ; 
            adenineComplement = adenineComplement + 1 ;
            complementaryStrand = [complementaryStrand, 'A'] ; 
        elseif n == 'G'
            guanineCount = guanineCount + 1 ; 
            cytosineComplement = cytosineComplement + 1 ; 
            complementaryStrand = [complementaryStrand, 'C'] ; 
        elseif n == 'C'
            cytosineCount = cytosineCount + 1 ; 
            guanineComplement = guanineComplement + 1 ; 
            complementaryStrand = [complementaryStrand, 'G'] ; 
        end 
    end 
end 

% (H) Similar logic and function to above, empty char with for-loop to iterate through and create mRNA strand
mRNAStrand = '' ; 
for i = 1:length(codingStrand) 
    n = codingStrand(i) ;
    if n == 'A'
        mRNAStrand = [mRNAStrand, 'A'] ;
    elseif n == 'T'
        mRNAStrand = [mRNAStrand, 'U'] ;
    elseif n == 'G'
        mRNAStrand = [mRNAStrand, 'G'] ;
    elseif n == 'C'
        mRNAStrand = [mRNAStrand, 'C'] ;
    end
end

% Learned about dictionaries in Python, MATLAB has containers, so I used this to map out the codon chart, so I made two arrays and each element in the arrays correspond 1:1 so codon is matched with respective amino acid then initialized the container as a map with codons as keys and amino acids as the value associated with them (Python dictionaries were less of a headache)
codonList =  ["UUU", "UUC", "UUA", "UUG", "UCU", "UCC", "UCA", "UCG", "UAU", "UAC", "UGU", "UGC", "UGG", "CUU", "CUC", "CUA", "CUG", "CCU", "CCC", "CCA", "CCG", "CAU", "CAC", "CAA", "CAG", "CGU", "CGC", "CGA", "CGG", "AUU", "AUC", "AUA", "ACU", "ACC", "ACA", "ACG", "AAU", "AAC", "AAA", "AAG", "AGU", "AGC", "AGA", "AGG", "GUU", "GUC", "GUA", "GUG", "GCU", "GCC", "GCA", "GCG", "GAU", "GAC", "GAA", "GAG", "GGU", "GGC", "GGA", "GGG", "AUG",  "UAA",  "UAG",  "UGA"] ;
aminoAcids = ["Phe", "Phe", "Leu", "Leu", "Ser", "Ser", "Ser", "Ser", "Tyr", "Tyr", "Cys", "Cys", "Trp", "Leu", "Leu", "Leu", "Leu", "Pro", "Pro", "Pro", "Pro", "His", "His", "Gln", "Gln", "Arg", "Arg", "Arg", "Arg", "Ile", "Ile", "Ile", "Thr", "Thr", "Thr", "Thr", "Asn", "Asn", "Lys", "Lys", "Ser", "Ser", "Arg", "Arg", "Val", "Val", "Val", "Val", "Ala", "Ala", "Ala", "Ala", "Asp", "Asp", "Glu", "Glu", "Gly", "Gly", "Gly", "Gly", "Met", "STP", "STP", "STP"] ; % "STOP" was abbreviated to "STP"
codonChart = containers.Map(codonList, aminoAcids) ; 

% Initialized an empty array to store codons which were defined by parsing and iterating the mRNA strand with a for-loop and looking at every three base (one codon) and grabbed each codon and store it in a "codon" and store the codon in the array called "codons"
codons = []; 
for i = 1:3:length(mRNAStrand)-2
    codon = mRNAStrand(i:i+2) ; 
    codons{end+1} = codon ; 
end 

% (I) Initialize empty char "aminoAcidSequence" and using a for-loop accessed the "codons" array and the value (sequence) of each codon, referenced it to check if it was in the "codonChart" Map and then set the amino acid from the "aminoAcids" array equal to aminoAcid and appended it to the aminoAcidSequence 
aminoAcidSequence = '' ;
for i = 1:length(codons)
    codon = codons{i} ;    
        if isKey(codonChart, codon)
        aminoAcid = codonChart(codon) ; 
        aminoAcidSequence = [aminoAcidSequence, aminoAcid] ; 
    end 
end

% (J) Okay stay with me for this one, again, initialize a container as a map, however this time I have a for-loop to create the container, the for-loop iterates through the "aminoAcidSequence" in sets of three (like previously with the mRNA) and checks if the amino acid is already a key in the container, if yes then we just increase the count + 1, if not, it creates the key and sets the count = 1 
aminoAcidCount = containers.Map('KeyType', 'char', 'ValueType', 'double') ; 
for i = 1:3:length(aminoAcidSequence)-2
    aminoAcid = aminoAcidSequence(i:i+2) ; 
    if isKey(aminoAcidCount, aminoAcid) 
        aminoAcidCount(aminoAcid) = aminoAcidCount(aminoAcid) + 1 ; 
    else
        aminoAcidCount(aminoAcid) = 1 ; 
    end 
end 

% (E) Just simple arithmetic to get counts and display/print information 
dsAdenine = adenineCount + adenineComplement ; 
dsThymine = thymineCount + thymineComplement ; 
dsGuanine = guanineCount + guanineComplement ; 
dsCytosine = cytosineCount + cytosineComplement ; 
nucleotideTotal = dsAdenine + dsThymine + dsGuanine + dsCytosine ; 

fprintf("\n (E) Total (Double Stranded) Adenines: %d\n", dsAdenine)
fprintf("\n (E) Total (Double Stranded) Thymines: %d\n", dsThymine)
fprintf("\n (E) Total (Double Stranded) Guanines: %d\n", dsGuanine)
fprintf("\n (E) Total (Double Stranded) Cytosines: %d\n", dsCytosine)

% (F) Hydrogen Bonds
hydrogenBonds = (adenineCount * 2) + (thymineCount * 2) + (guanineCount * 3) + (cytosineCount * 3) ;
fprintf("\n (F) Number of Hydrogen Bonds: %d\n", hydrogenBonds) ;

% (G) Molecular Weight information and formulas 
mwAdenine = 313.21; 
mwThymine = 304.2; 
mwGuanine = 329.21; 
mwCytosine = 289.18; 
% units in g/mol
molecularWeight1 = (mwAdenine * dsAdenine) + (mwThymine * dsThymine) + (mwGuanine * dsGuanine) + (mwCytosine * dsCytosine) + 79.0 ; 
molecularWeight2 = (nucleotideTotal * 607.4) + 157.9 ;
fprintf("\n (G) Molecular Weight of this particular dsDNA: %d", molecularWeight1)
fprintf(" g/mol or Da \n")

fprintf("\n     Input Strand Type: %s\n", strandType)
fprintf("\n     Inputted Strand: %s\n", dnatext)

% (C) Printing Complementary Strand
fprintf("\n (C) Complementary Strand: %s\n", complementaryStrand)

% I kept track of template and coding strand just so I can double check my ouputs easier
fprintf("\n     Template Strand (3' to 5'): %s\n", templateStrand)
fprintf("\n     Coding Strand (5' to 3'): %s\n", codingStrand)

% (H) Print out Amino Acid Sequence 
fprintf("\n (H) mRNA Strand (5' to 3'): %s\n", mRNAStrand)
% (I) Primary Amino Acid Chain/Strand
fprintf("\n (I) Amino Acid Sequence: %s\n", aminoAcidSequence)

% (J) Count and print out total occurance of each Amino Acid, for-loop goes through the keys in "aminoAcidCount" which are the amino acids with their corresponding count value, then print the amino acid followed by their count, loop iterates for all amino acids in the container 
fprintf("\n (J) Amino Acid Counts: \n")
keys = aminoAcidCount.keys ; 
for i = 1:length(keys)
    aminoAcid = keys{i} ; 
    count = aminoAcidCount(aminoAcid) ; 
    fprintf("     %s: %d\n", aminoAcid, count) ;
end

% Data presentation, use subplot so I can display both plots 
subplot(2,1,2) ; 
data = [ adenineCount adenineComplement dsAdenine ; thymineCount thymineComplement dsThymine ; guanineCount guanineComplement dsGuanine; cytosineCount cytosineComplement dsCytosine] ; 
bar(data)
ylabel("Amount of Base")
xlabels = ["Adenine" ; "Thymine" ; "Guanine" ; "Cytosine"] ; 
legend("Inputted Strand","Complementary Strand", "Double Strand")
set(gca, "XTickLabel", xlabels)
title(" (D) Base Counts") ; 

% Gotta convert the cell array to a matrix because the bar() function cannot take in a cell array to plot
subplot(2,1,1) ; 
values = cell2mat(aminoAcidCount.values) ; 
bar(values)
xticks(1:length(keys)) ; 
xticklabels(keys) ; 
ylabel("Amino Acid Count")
xlabel("Amino Acid")
title(" (J) Amino Acid Counts")
sgtitle(" Amino Acid and Base Counts")

% (N) References: 
% MW Formula 1 : https://www.thermofisher.com/us/en/home/references/ambion-tech-support/rna-tools-and-calculators/dna-and-rna-molecular-weights-and-conversions.html
% MW Formula 2 : https://www.sciencedirect.com/topics/nursing-and-health-professions/molecular-weight#:~:text=Anhydrous%20molecular%20weight%20of%20each,salt%20solution)%20is%20325%20Daltons.
% General MATLAB Help: 
% (How to use containers.Map) 
% https://stackoverflow.com/questions/56804013/create-a-matlab-dictionary-like-in-python (I took Python, Python > MATLAB)
% https://www.mathworks.com/help/matlab/ref/containers.map.html
% https://www.mathworks.com/help/matlab/ref/containers.map.keys.html
% https://www.mathworks.com/help/matlab/ref/cell2mat.html

% Collaborated with : Erick Garcia and Nohemi Rodriguez Hernandez