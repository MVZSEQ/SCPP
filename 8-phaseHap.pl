#phase haplotypes
#need to install: GATK, picard-tools, SamTools
#version 1.09 June 28 2013
#Ke Bi-> kebi@berkeley.edu

use strict; 
use warnings;
use Getopt::Std;
use File::Basename;


&main;
exit;

sub main {
        &usage if (@ARGV<1);
        my $command = shift(@ARGV);
        my %fun = (phasing=>\&phasing);
        die("Unknown command \"$command\"\n") if (!defined($fun{$command}));
        &{$fun{$command}};
      }

sub usage {
  die(qq/
Usage:  phaseHap.pl <command> [<arguments>]\n
Command: phasing:  Haplotype phasing 
\n/);
}


  
sub phasing {
  die(qq/
Usage phaseHap.pl phasing [options]

Options: -G CHR  GATK folder
         -p CHR  Picard folder
         -b CHR  Bam files folder
         -f CHR  Reference files folder
	 -r CHR  Where results will go
	 -D INT  Minimal coverage filter [10]
	 -Q INT  Minimal quality filter [21]
	 -M FLOAT Minimal MAF filter [0.2]
         -H FLOAT Minimal freq to call a homozygous site [0.8]

Note: folders do not end with "\/"; 

\n/) if (!@ARGV);
  
  
  my %opts = (G=>undef, p=>undef, b=>undef,f=>undef, r=>undef, D=>10, Q=>21, M=>0.2, H=>0.8);
  getopts('G:p:b:f:r:D:Q:M:H:', \%opts);
  
  my $ref_dir = $opts{f}. "/" ;  
  my @bams = <$opts{b}/*sorted.bam>;  
  my $gatk = $opts{G}.'/GenomeAnalysisTK.jar';
  my $resultsDir = $opts{r}. "/";
  my $AddOrReplace = $opts{p}.'/AddOrReplaceReadGroups.jar'; 
  my $CreateSequenceDictionary = $opts{p}.'/CreateSequenceDictionary.jar';
   
  foreach my $bams  (@bams) {    
    my $ref = $ref_dir. $1 .'.fa' if basename($bams) =~ m/(\S+).sorted.bam/;
    system ("samtools faidx $ref");	
    my $lib = $1  if basename($bams) =~ m/(\S+).sorted.bam/;	
    my $index_ref = $ref_dir. $lib . '.dict';    
    system("java -jar $AddOrReplace INPUT=$bams OUTPUT=$resultsDir$lib.rg.bam RGID=$lib RGLB=exon_capture RGPL=illumina RGPU=lane3 RGSM=$lib");   
    system ("java -jar $CreateSequenceDictionary R=$ref O=$index_ref");    
    system("samtools sort $resultsDir$lib.rg.bam $resultsDir$lib.rg.sort");    
    system("samtools index $resultsDir$lib.rg.sort.bam");    
    system("java -Xmx8g -jar $gatk -T HaplotypeCaller -R $ref -I $resultsDir$lib.rg.sort.bam -o $resultsDir$lib.rg.vcf");    
    system("java -Xmx8g -jar $gatk -T ReadBackedPhasing -R $ref -I $resultsDir$lib.rg.sort.bam --variant $resultsDir$lib.rg.vcf --min_base_quality_score $opts{Q} -o $resultsDir$lib.raw.vcf");
    system("java -Xmx8g -jar $gatk -T VariantAnnotator -A DepthPerAlleleBySample -A AlleleBalance -A FisherStrand -A HaplotypeScore -A HardyWeinberg -R $ref -I $resultsDir$lib.rg.sort.bam --variant $resultsDir$lib.raw.vcf -o $resultsDir$lib.annotated.vcf");
    system ("grep -v \'#\'  $resultsDir$lib.annotated.vcf > $resultsDir$lib.raw.vcf_org");
    system ("grep -v \'ABHet=1.00\' $resultsDir$lib.raw.vcf_org > $resultsDir$lib.raw2.vcf ");
    system ("rm  $resultsDir$lib.raw.vcf* $resultsDir$lib.annotated.vcf* $resultsDir$lib.rg.bam* $resultsDir$lib.rg.sort* $resultsDir$lib.rg.vcf* ");
    open (RAW, "<", $resultsDir . $lib.'.raw2.vcf');
    system ("rm $resultsDir$lib.raw2.vcf* ");
    
    chomp (my @first_line = split /\s/, <RAW>); 
    
    open (FILE, ">", $resultsDir .  $first_line[0]. '.vcf' ); 
    print FILE join ("\t", @first_line), "\n";
    
    foreach (<RAW>) {
      
      chomp (my @line = split /\t+/, $_);
      if (($line[0] !~ m/^#/)) {
	if ($line[0] eq $first_line[0]) {
	  print FILE join ("\t", @line), "\n";
	  
	}
	if ($line[0] ne $first_line[0]) {
          open (FILE, ">", $resultsDir . $line[0]. '.vcf' );
	  print FILE join ("\t", @line), "\n";
	  @first_line  =   @line;
	  
	}
      }
    }
    close FILE;
    close RAW;
    
    
    
    my @new_files = <$resultsDir*vcf>;
    
    foreach my $file (@new_files) {
      
      my $a = 0;
      
      my $contig = $1 if (basename($file) =~ m/(\S+).vcf/);
      
      open (NUM, ">", $resultsDir . $contig. '.line'); 
      
      open (EACH, "<", $file);
      
      
      foreach (<EACH>) {
	my @line = split /\t+/, $_;	
	if (($line[0] !~ m/^#/) && (substr($line[9],0,3) eq '0/1') ) {
	  print NUM $line[1], "\n";
	  $a ++;
	}
      }
      
      close NUM;
      close EACH;
      
      
      if ($a >1) {
	
	open (TINA, "<", $resultsDir. $contig . '.line');
	
	
	my @array1;
	
	foreach (<TINA>) {
	  chomp($_);
	  push (@array1, $_);  
	}
	my $size = scalar (@array1);
	close TINA;
	
	
	for (my $b = 1; $b <  $size; $b++) {
	  
	  
	  
	  my $newraw = $resultsDir . $contig . '.vcf'. '_' .  $b ;
	  open (NEWRAW, ">", $newraw);
	  
	  my $raw_final = $resultsDir . $contig . '.vcf' . '_' . $size ;
	  open (RAW_final,">", $raw_final);
	  
	  my $raw = $resultsDir . $contig. '.vcf';
	  open (RAW1, "<", $raw);
	  
	  foreach (<RAW1>) {
	    chomp(my @line = split /\t+/, $_);
	    if ($b == 1 && $b < $size-1) { 
	      if (($line[1] < $array1[$b])) {
		print NEWRAW join ("\t", @line), "\n";
	      }
	    }
	    
	    
	    if  ($b == 1 && $b == $size-1) { 
	      if (($line[1] < $array1[$b])) {
		print NEWRAW join ("\t", @line), "\n";
	      }
	      if (($line[1] >= $array1[$b])) {
		print RAW_final join ("\t", @line), "\n";
	      } 
	    }
	    
	    
	    if ($b > 1 && $b < $size-1) {
	      if (($line[1] < $array1[$b]) && ($line[1] >= $array1[$b-1]) ) {
		print NEWRAW join ("\t", @line), "\n";	      
	      }
	      
	    }
	    
	    if ($b > 1 && $b == $size-1) {
	      if (($line[1] < $array1[$b]) && ($line[1] >= $array1[$b-1]) ) {
		print NEWRAW join ("\t", @line), "\n";	      
	      }
	      if (($line[1] >= $array1[$b])) {
		print RAW_final join ("\t", @line), "\n";
	      } 
	    }
	    
	  } #foreach (<RAW>)
	  
          close NEWRAW;
	  close RAW_final;   
          close RAW1;
	  
	} #for (my $b = 0)
	system ("rm $resultsDir$contig.vcf");
      }  #if a>1    
      
      
      if ($a<=1) {
	system ("mv $resultsDir$contig.vcf $resultsDir$contig.vcf_0");
      }
      
    } # foreach my $file (@new_files)    
    
    
    my %ref;
    open (REF, "<", $ref);
    while (<REF>) {
      chomp (my $line = $_);
      if ($line =~ m/^>([A-Z|a-z|0-9]+)/) {
	my $ID = $1;
	chomp(my $seq = <REF>);
        $ref{$ID} = $seq;
      }
    }
    
    close REF;
    
    my @final_line = <$resultsDir*line>;
       
    my $new_ref = $ref_dir . $lib . '_new.fa';
    open (NEWREF, ">",$new_ref);
    foreach my $file (@final_line) {
      open (ANNA, "<", $file);
      my @array2;
      foreach (<ANNA>) {
	chomp($_);
	push (@array2, $_);  
      }
      my $size = scalar (@array2);
      close ANNA;
      
      if ($size <= 1) {
	
	my $contig = $1 if (basename($file) =~ m/(\S+).line/);
	print NEWREF ">", $contig, "\n", $ref{$contig}, "\n"; 
      }
      elsif ($size > 1) {
	my $contig = $1 if (basename($file) =~ m/(\S+).line/);
	open (JES, "<", $file); 
	chomp (my $first = <JES>);
	my $new_size = $size - 1; #=2
	
	my $a = 0;
	my $b = 1;
	my $totlen = length ($ref{$contig});
	
	foreach (<JES>) {
	  chomp (my @line = split /\t/, $_);
	  
	  if ($b < $new_size){
	    print NEWREF ">", $contig, "_vcf_", $b, "\n",  'N' x $a, substr ($ref{$contig}, $a , $line[0]-$a-1), 'N' x ($totlen-$a-($line[0]-$a-1)), "\n";
	    open (VCF, "<",  $resultsDir . $contig . '.vcf_' . $b);
	    open (NEWVCF, ">", $resultsDir . $contig . '.vcf2_' . $b);
	    foreach (<VCF>) {
	      chomp (my @l = split /\t+/, $_);
	      print NEWVCF $contig, "_vcf_", $b, "\t", join ("\t", @l[1 .. $#l]), "\n";
	    }
	    close VCF; close NEWVCF;
	    
	    system ("rm  $resultsDir$contig'.vcf_'$b ");
	    $a = $line[0]-1;
	    $b++;
	  }
	  
	  elsif ($b = $new_size) {
	    print NEWREF ">", $contig, "_vcf_",$b, "\n", 'N' x $a, substr ($ref{$contig}, $a , $line[0]-$a-1), 'N' x ($totlen-$a-($line[0]-$a-1)),  "\n";
	    open (VCF, "<", $resultsDir . $contig . '.vcf_' . $b);
	    open (NEWVCF, ">", $resultsDir . $contig . '.vcf2_' . $b );
	    foreach (<VCF>) {
	      chomp (my @l = split /\t+/, $_);
	      print NEWVCF $contig, "_vcf_", $b, "\t", join ("\t", @l[1 .. $#l]) , "\n" ;
	    }
	    close VCF; close NEWVCF;
	    system ("rm  $resultsDir$contig'.vcf_'$b ");
	    
	    $a = $line[0]-1; 
	    print NEWREF ">", $contig, "_vcf_", $b+1, "\n",'N' x $a,  substr ($ref{$contig}, $a), "\n";
	    my $d = $b+1;
	    open (VCF, "<", $resultsDir . $contig . '.vcf_' . $d);
	    open (NEWVCF, ">", $resultsDir . $contig . '.vcf2_' . $d );
	    foreach (<VCF>) {
	      chomp (my @l = split /\t+/, $_);
	      print NEWVCF $contig, "_vcf_", $d, "\t", join ("\t", @l[1 .. $#l]) , "\n";
	    }
	    system ("rm  $resultsDir$contig'.vcf_'$d ");
	    close VCF; close NEWVCF;
	  } 
	}
      }
    }
    close NEWREF;
    
    system ("cat $resultsDir*vcf* > $resultsDir$lib'.final' ");
    system ("rm $resultsDir*line $resultsDir*vcf*");
    
    
    
    
    open (REF, "<", $new_ref);
    my %seq; my %seqa; my %seqb;
 
     while (<REF>) {
      chomp();
      my $gene = $1 if ($_ = m/^>(\S+)/);
      chomp (my $second = <REF>);
      $seq{$gene} = $second;
      $seqa{$gene} = $second;
      $seqb{$gene} = $second;
    }
    close REF;
    
    
    
    
    my $hap1 = $resultsDir.$lib.'.hap1.txt';
    my $hap2 = $resultsDir.$lib.'.hap2.txt';
    
    open (HAP1, ">", $hap1);
    open (HAP2, ">", $hap2);
    
    open (FINAL, "<", $resultsDir. $lib . ".final");
    open (OUTVCF, ">", $resultsDir. $lib . ".final2");
    
    my @VCF = reverse <FINAL>;
    
    foreach (@VCF) {
      chomp (my @line = split /\t+/, $_);
      print OUTVCF join ("\t", @line), "\n";
    }
    close OUTVCF;
    close FINAL;
    
    system ("rm $resultsDir$lib'.final' ");
    
    open (NEW, "<",  $resultsDir. $lib . ".final2");
    chomp(my @first1 = split /\t+/, <NEW>);
    
    close NEW;
    
    my $first = $first1[0];
    my $seq1 ;
    my $seq2 ;
    
    open (NEW2, "<", $resultsDir. $lib . ".final2" );
    system ("rm $resultsDir$lib'.final2' ");
    while  (<NEW2>) {
      chomp(my @line = split /\t/, $_);
      my $gene = $line[0];
      my $pos = $line[1];
      my $depth = $1 if ($line[7] =~ m/;DP=(\d+);/); 
      my $refbase = $line[3];
      my $altbase = $line[4];
      my $state = substr($line[7],0,5) if ($line[7] =~ m/^ABH/) ;
      my $freq = substr($line[7],6,4) if ($line[7] =~ m/^ABH/) ;
      my $call_all = substr($line[9],0,3);
      my $len_ref = length ($refbase);
      my $len_alt = length ($altbase);
      #my $ref_count = $1 if ($line[9] =~ m/\d[\/|\\|\|]\d:(\d+),(\d+):\d+:/ );
      #my $alt_count = $2;
      #my $indel_depth = $ref_count + $alt_count; 
           
      my $a1 = 'N';
      my $a2 = 'N';    
     
      
      if ($seq{$gene}) { 	
	if ($depth >= $opts{D}) { 		  
	  if ($line[7] =~ m/^ABH/) {
	    if ($state eq 'ABHet') {	      
	      if ($freq >= $opts{M}) {  
		if (($state eq 'ABHet') && ($call_all eq '0/1')) {
		  $a1 = lc ($refbase); 
		  $a2 = lc ($altbase);
		}
		if (($state eq 'ABHet') && ($call_all eq '0|1')) {
		  $a1 = $refbase; 
		  $a2 = $altbase;	    
		}
		if (($state eq 'ABHet') && ($call_all eq '1|0')) {
		  $a1 = $altbase; 
		  $a2 = $refbase;
		}	
	      }
	      
	      if ($freq < $opts{M}) {
		$a1 = $a2 = 'N';
	      }
	    }
	    
	    if ($state eq 'ABHom') { 
	      if (($freq >= $opts{H}) && ($call_all eq '1|1')) {
		$a1 = $a2 = $altbase;
	      }
	      if (($freq < $opts{H}) && ($call_all eq '1|1')) {
		$a1 = $a2 = 'N';
	      }
	    }    
	  } #if ($line[7] =~ m/^ABH/ )
	  
	  
	  if ($line[7] !~ m/^ABH/) {
	    if ($altbase =~ m/,/) {
	      $a1 = $a2 = 'N';
	    } 	    
	    if ($altbase !~ m/,/) {
	      if (($call_all eq '0/1') || ($call_all eq '0|1')) {
		$a1 = $refbase;
		$a2 = $altbase;
	      }
	      if (($call_all eq '1/0') || ($call_all eq '1|0')) {
		$a1 = $altbase;
		$a2 = $refbase;
	      }
	      if (($call_all eq '1|1') || ($call_all eq '1/1')) {
		$a1 = $a2 = $altbase;
	      }  
	    } 
	  }
	} #($depth >= $opts{D}) 
	
	
	
	if ($depth < $opts{D}) { 
	  if ($line[7] =~ m/^ABH/) {
	    $a1 = $a2 = 'N';
	  }
	  if ($line[7] !~ m/^ABH/) { 
	    if ($altbase =~ m/,/) { 
	      $a1 = $a2 = 'N' x $len_ref;	      
	    }
	    
	    if ($altbase !~ m/,/) { 
	      if (($call_all eq '0/1') || ($call_all eq '0|1')) {
		$a1 = 'n' x $len_ref;
		$a2 = 'n' x $len_alt;
	      }
	      if (($call_all eq '1/0') || ($call_all eq '1|0')) {
		$a1 = 'n' x $len_alt;
		$a2 = 'n' x $len_ref;
	      }
	      if (($call_all eq '1|1') || ($call_all eq '1/1')) {
		$a1 = $a2 = 'n' x $len_alt;
	      } 
	    }
	  }
	} # ($depth < $opts{D})
	
      } # if ($seq{$gene})
      
      
      
      if ($gene eq $first) {
	if (!eof(NEW2) ) {	
	  $seq1 = $seqa{$first};
	  $seq2 = $seqb{$first};
	  substr($seq1, $pos-1, $len_ref) = $a1;
	  substr($seq2, $pos-1, $len_ref) = $a2;
	  $seqa{$first} = $seq1;	  
	  $seqb{$first} = $seq2;	  
	  $first = $gene;
	}
	if (eof (NEW2)) {	
	  $seq1 = $seqa{$first};
	  $seq2 = $seqb{$first};
	  substr($seq1, $pos-1, $len_ref) = $a1;
	  substr($seq2, $pos-1, $len_ref) = $a2;
	  print HAP1 ">", $first, "\n", $seq1, "\n";
	  print HAP2 ">", $first, "\n", $seq2, "\n";
	}
      }
      
      if ( $gene ne $first ) {
	if (!eof(NEW2) ) { 
	  print HAP1 ">", $first, "\n", $seq1, "\n";
	  print HAP2 ">", $first, "\n", $seq2, "\n";	
	  $seq1 = $seqa{$gene};
	  $seq2 = $seqb{$gene};
	  substr($seq1, $pos-1, $len_ref) = $a1;
	  substr($seq2, $pos-1, $len_ref) = $a2;
	  $first = $gene;
	  $seqa{$first} = $seq1;	  
	  $seqb{$first} = $seq2;
	  
	}
      	if (eof (NEW2)) {	
	  print HAP1 ">", $first, "\n", $seq1, "\n";
	  print HAP2 ">", $first, "\n", $seq2, "\n";	
	  $seq1 = $seqa{$gene};
	  $seq2 = $seqb{$gene};
	  substr($seq1, $pos-1, $len_ref) = $a1;
	  substr($seq2, $pos-1, $len_ref) = $a2;
	  print HAP1 ">", $gene, "\n", $seq1, "\n";
	  print HAP2 ">", $gene, "\n", $seq2, "\n";
	}
	
      }
      
      
      
    }
  
    close HAP1;
    close HAP2;
    
    my $homo = $resultsDir.$lib.'.homo_hap.txt';
    open (HOMO, ">", $homo);
    foreach my $d (keys %seq) {
      if ($d =~ m/(\S+)_vcf_\d+/) {
	delete $ref{$1} if $ref{$1};       
      }
      else {
	delete $ref{$d} if $ref{$d};        
      }
    }
  
  
  
  foreach my $c (keys %ref) {
    print HOMO ">", $c, "\n", $ref{$c}, "\n";
  }
  close HOMO;
  
  my $final_hap1 = $resultsDir.$lib.'.final_hap1.txt';
  my $final_hap2 = $resultsDir.$lib.'.final_hap2.txt';
    
    system ("cat $hap1 $homo >  $final_hap1");
    system ("cat $hap2 $homo >  $final_hap2");
    unlink ($hap1,$hap2,$homo);
    system ("rm $ref_dir$lib'.dict'  $ref_dir$lib'.fa.fai'");
    
    
    
    
  } # foreach my $bams  (@bams) 
} #sub phasing
