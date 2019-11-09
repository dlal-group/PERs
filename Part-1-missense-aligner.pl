#!/usr/bin/env perl
#use warnings;
use strict;
use Data::Dumper;
use List::MoreUtils qw/ uniq /;

########################################################### 
print "STEP 1\n\tIndexing alignment\n\tAnnotating paralog score\n\tCreating input files\n";
###########################################################

# ====================== Hash parazscore
open PARAZ, "<./db/hg19.paralog.allScores" or die "$!";

my @para_array = <PARAZ>;
my %parazscore_hash;
foreach (@para_array) {
	chomp($_);
	my ($key, $parazscore) = split /\t/, $_;
	$parazscore_hash{$key} = $parazscore;
}
close PARAZ;

########################################################### Loading and formating MUSCLE output
my @files_clw = glob '*.clw' or die "$!";
foreach my $archivo1 (@files_clw) {
	chomp($archivo1);
	open RAW,"<$archivo1" or die "$!";
	$archivo1 =~ s/\.clw//gi;
	open RESULTADOS1, ">$archivo1.raw.out" or die "$!";
	open CURATED, ">$archivo1.clw.curated" or die "$!";	
	my @raw = <RAW>;
	my @curated = ();
	shift(@raw);
	shift(@raw);
	foreach (@raw) {
		$_ =~ s/ +/ /g;
		$_ =~ s/\A\n//g;
		$_ =~ s/:|\.|\*//g;
		$_ =~ s/ +\n//g;
		$_ =~ s/ /\t/g;
		push (@curated, $_);
	}
	my @col1 = ();
	foreach (@curated) {
		my $gene = (split /\t/, $_)[0];
		push (@col1, $gene);
		print CURATED "$_";
	}
	close CURATED;

########################################################### Transpose alignment and indexing
	
 	open ALIGN, "<$archivo1.clw.curated" or die "$!";
 	my @align = <ALIGN>;
 	my @genes = uniq @col1;
 	my $members = $#genes;
 	my $loops = (($#align+1)/$members);
 	my @align_lineal = ();
 	# ====================== Making continuous sequence 
 	for my $m (0..($members-1)) {
 		my $seq_cat = ();
 		my $seq = ();
 		my $gene = (split /\t/, $align[$m])[0];
 		for my $c (1..$loops) {
 			$seq = (split /\t/, $align[$m])[1];             
 			chomp($seq);
 			$seq_cat = $seq_cat.$seq;
 			$m = $m + $members;
 			if ($c eq $loops)  {
 				push @align_lineal, $gene."\t".$seq_cat;
 			}
 		}
 	}
 	# ================ Split in tabs 
 	my @align_tab = ();                                    
 	foreach (@align_lineal) {
 		my ($gid, $lseq) = split /\t/, $_;
 		my $tabseq = join("\t", split(//, $lseq));          
 		push @align_tab, $gid."\t".$tabseq;
 	}
 	# ================ Add Gene_Index
 	my @align_num = ();
 	foreach (@align_tab) {
 		my @linea_num = ();
 		my @linea = split /\t/, $_;
 		my $gi = 0; 
 		for my $pos (0..$#linea) {
 			if ($pos == 0) {
 				$linea_num[$pos] = $linea[$pos];
 			} else {
 				if ($linea[$pos] eq "-") {
 					$linea_num[$pos] = $linea[$pos];					
 				} else {
 					$gi++;
 					$linea_num[$pos] = $linea[$pos]."_".$gi;
 				}
 			}
 		}
 		my $linea_out = join "\t", @linea_num;
 		unshift @align_num, $linea_out;
 	}
 	# ================ Add Family Index        
 	my @toindex = split /\t/, $align_num[0];
 	my @conteo = (0..$#toindex);
 	$conteo[0] = "Index";
 	my $index = join "\t", @conteo;
 	unshift @align_num, $index;
 	# ================ Transpose                         
 	my @transpuesta = ();
 	for my $i (0..$#align_num) {                            
 		my $fila = shift @align_num;                        
 		chomp($fila);                                       
 		my @columnas = split /\t/, $fila;                   
 		my $ref_columnas = \@columnas;                      
 		foreach my $c (0..$#columnas) {                     
 			push @{$transpuesta[$c]}, $ref_columnas->[$c];  
 		}                                                   
 	}
 	# ================ Print	
 	for my $nueva_fila (@transpuesta) {                     
 	  for my $nueva_col (@{$nueva_fila}) {                  
 	      print RESULTADOS1 $nueva_col, "\t";               
 	  }
 	  print RESULTADOS1 "\n";                               
 	}

 	# ================ Closing statements
 	close RAW;
 	close ALIGN;
 	close RESULTADOS1;
 	open RESULTADOS2, "<$archivo1.raw.out" or die "$!";
 	open RESULTADOS3, ">$archivo1.out" or die "$!";
 	my @res = <RESULTADOS2>;	
 	foreach (@res) {
 		$_ =~ s/\t\n/\n/g; 
 		print RESULTADOS3 "$_";
 	}
 	close RESULTADOS2;
 	close RESULTADOS3;
 	unlink("$archivo1.clw.curated") or die "Could not delete the file!\n";
 	unlink("$archivo1.raw.out") or die "Could not delete the file!\n";

############################################################################### Add parazscore
  	
 	open RESULTADOS4, "<$archivo1.out" or die "$!";
 	my @input = <RESULTADOS4>;
 	my $header = shift(@input);
 	chomp($header);
 	my @genes_in_order = split /\t/, $header;
  	open OUT, ">$archivo1.input";
 	my @output = ();
 	unshift(@output, $header."\tParazscore");
 	foreach (@input) {
 		my $fila_out = $_;
 		chomp($fila_out);
 		my @fila = split /\t/, $fila_out;
 		my $exito = 0;
 		for my $g (1..$#fila) {                        #
 			my $keyp = $genes_in_order[$g]."_".$fila[$g];
 			chomp($keyp);
 			if (exists $parazscore_hash{$keyp}) {
 				$fila_out = $fila_out."\t".$parazscore_hash{$keyp};
 				$exito = 1;
 				last;
 			}
 		}
 		if ($exito == 1) {
 			push @output, $fila_out;	
 		} else {
 			push @output, $fila_out."\tNA";
 		}
	
 	}
 	foreach (@output) {
 		print OUT $_."\n";
 	}
 	close OUT;
 	unlink("$archivo1.out") or die "Could not delete the file!\n";
}

################################################################## Working with gene-only fasta files
my @files_fasta = glob '*.fasta' or die "$!";
foreach my $fasta (@files_fasta) {
	chomp($fasta);
	open FASTA, "<$fasta" or die "$!";
	$fasta =~ s/\.fasta//gi;
	open GENEWISE, ">$fasta.input" or die "$!"; 
	my @genefasta = <FASTA>;
	my $genename = shift(@genefasta);
	$genename =~ s/\A>//gi;
	chomp($genename);
	my $fastaseq = join //,@genefasta;
	$fastaseq =~ s/\n//gi;
	chomp($fastaseq);
	print GENEWISE "Index\t$genename\tParazscore\n";
	my @aa = split //, $fastaseq;
	my $cont = 0;
	foreach my $a (@aa) {
		chomp($a);
		$cont = $cont +1;
		my $query = $genename."_".$a."_".$cont;
		if (exists $parazscore_hash{$query}) {
			print GENEWISE $cont."\t".$a."_".$cont."\t".$parazscore_hash{$query}."\n";
		} else {
			print GENEWISE $cont."\t".$a."_".$cont."\tN/A\n";
		}
	}
	close GENEWISE;
	close FASTA;
}

########################################################### 
print "STEP 2\n";
###########################################################

# ====================== Genetic code hash
open ARCHIVO1, "<./db/aminoacid-code" or die "$!";
my @array = <ARCHIVO1>;
my %cualaa;
for my $h (0..$#array) {
	my $fila = shift @array;
	chomp($fila);
	my ($key, $letra) = (split /\t/, $fila)[2,0];
	chomp($key, $letra);
	$cualaa{$key} = $letra;
}
my %revcualaa = reverse %cualaa;
my $n_mutaciones = 0;

# ====================== Loading clinvar and gnomAD missense variants
print "\tLoading input files\n";
open ARCHIVO2, "<./db/input.clinvar-hgmd" or die "$!";
open ARCHIVOG1, "<./db/input.gnomad.1" or die "$!";
open ARCHIVOG2, "<./db/input.gnomad.2" or die "$!";
open ARCHIVOG3, "<./db/input.gnomad.3" or die "$!";
open ARCHIVOG4, "<./db/input.gnomad.4" or die "$!";
open ARCHIVOG5, "<./db/input.gnomad.5" or die "$!";
open ARCHIVOG6, "<./db/input.gnomad.6" or die "$!";
open ARCHIVOG7, "<./db/input.gnomad.7" or die "$!";
my @clinvarhgmd = <ARCHIVO2>;
my @gnomad1 = <ARCHIVOG1>;
my @gnomad2 = <ARCHIVOG2>;
my @gnomad3 = <ARCHIVOG3>;
my @gnomad4 = <ARCHIVOG4>;
my @gnomad5 = <ARCHIVOG5>;
my @gnomad6 = <ARCHIVOG6>;
my @gnomad7 = <ARCHIVOG7>;
my @gnomadall = ();
push(@gnomadall,@gnomad1);
push(@gnomadall,@gnomad2);
push(@gnomadall,@gnomad3);
push(@gnomadall,@gnomad4);
push(@gnomadall,@gnomad5);
push(@gnomadall,@gnomad6);
push(@gnomadall,@gnomad7);

# ====================== Main loop

for my $s (0..1) {
	my $loop = 0;
	# Load all input files: output os step 1
	my @files_input = glob '*.input' or die "$!"; 												#To GitHub
	my @mut_source = ();
	my $source = ();
	# Select variation source
	if ($s==1) {
		@mut_source = @gnomadall;
		$source = "gnomad";
	} else {
		@mut_source = @clinvarhgmd;
		$source = "clinvar-hgmd";
	}
	print "\tAnnotating $source missense variants\n";
	# Loading single input file 
	foreach my $archivo1 (@files_input) {
		chomp($archivo1);
		open ARCHIVO3,"<$archivo1" or die "$!";
		$archivo1 =~ s/\.input//;
		open OUT, ">./db/$archivo1.$source.binary";
		my @alineamiento = <ARCHIVO3>;
		my $header = $alineamiento[0];
		chomp($header);
		# Hash of all family members
		my %cualgen;
		my @genes = split /\t/, shift(@alineamiento);
		my $cont1 = 0;
		for my $g (0..$#genes) {
			my $nombre = $genes[$g];
			chomp($nombre);
			$cualgen{$nombre} = $cont1;
			$cont1=$cont1+1;
			$header= "$header\t$nombre";
		}
		$header =~ s/\tIndex//g;
		$header =~ s/\tParazscore\Z//g;
		# Creating empty final matrix "@alineamiento_final" full of 0s
		my @alineamiento_final = ();
		my $totalgenes = (scalar keys %cualgen)-2;                                            
		foreach my $fila (@alineamiento) {
			chomp($fila);
			my @elementos = split /\t/, $fila;
			for my $r (1..($totalgenes)) {                                                    
				if ($source =~ /gnomad/ ) {	
					if ($elementos[$r] =~ /-/) {
						push(@elementos,1);
					} else {
						push(@elementos,0);
					}
				} else  {
					push(@elementos,0);
				}
			}
			my $fila_final = join "\t", @elementos;
			push (@alineamiento_final, "$fila_final\t0\n");
		}
		# Nesting empty final matrix for future looping
		my @alineamiento_anidado;
		for my $f (0..$#alineamiento_final) {
			my $row = $alineamiento_final[$f];
			chomp($row);
			my @array2 = (split /\t/, $row);
			my $ref = \@array2;
			$alineamiento_anidado[$f] = $ref;
		}
		# Nesting disease elements
		my @alineamiento_disease;
		for my $f2 (0..$#alineamiento_final) {
			my $row2 = $alineamiento_final[$f2];
			chomp($row2);
			my @array3 = (split /\t/, $row2);
			my $ref2 = \@array3;
			$alineamiento_disease[$f2] = $ref2;
		}
		# Mapping missense variants from source

		for my $g (1..($#genes-1)) {
			chomp($genes[$g]);
			my @gene_mut_source = grep /\A($genes[$g])\t/, @mut_source;
			foreach my $mut (@gene_mut_source) {                              
				chomp($mut);                                                  
				my ($gen,$alelos,$aa3,$pos_mut,$disease) = (split /\t/, $mut)[0,3,7,9,10];
				chomp($disease);
				my $aa_ref = $cualaa{$aa3};
				my $col = $cualgen{$gen};
				my $cont2 = 0;
				my $pos_index = 0;
				my $aa_index_with_pos = ();
				my $aa_index = ();
				for my $i (0..$#alineamiento_final) {
					$aa_index_with_pos = $alineamiento_anidado[$i]->[$col];
					$aa_index = (split /_/, $aa_index_with_pos)[0];
					if (exists $revcualaa{$aa_index}) {
						$cont2 = $cont2 + 1;
						if ($cont2 eq $pos_mut) {
							$pos_index = $i;
							if ($aa_index eq $aa_ref) {
								$alineamiento_anidado[$i]->[$totalgenes+($col+1)] = "1";
								$alineamiento_disease[$i]->[$totalgenes+($col+1)] = "$gen:$disease";
								$n_mutaciones = $n_mutaciones +1;
							} else {
								# place to save variants that do not match 
							}
						}
					} 
				}
			}
		}
		# header of final matrix
		$header =~ s/\tIndex_[0-9]+\t/\t/;
		chomp($header); 
		$header = $header."\tTotal_Index\tGene:Disease\n";
		print OUT $header;
		
		# Annotating final matrix
		for my $k (0..$#alineamiento_final) {          
			for my $c (($totalgenes+2)..($totalgenes*2+1)) {
				if ($alineamiento_anidado[$k]->[$c] ne "-" ) {      # Porque no me cuenta el -? ver los outputs
					$alineamiento_anidado[$k]->[($totalgenes*2)+2] = $alineamiento_anidado[$k]->[($totalgenes*2)+2] + $alineamiento_anidado[$k]->[$c];
					$alineamiento_disease[$k]->[($totalgenes*2)+2] = $alineamiento_disease[$k]->[($totalgenes*2)+2].";".$alineamiento_disease[$k]->[$c];
				}
			}
		}
		# Printing final matrix to output file
		for my $y (0..$#alineamiento_final) {           
			my $escribir = join "\t", @{$alineamiento_anidado[$y]};
			my $all_diseases = $alineamiento_disease[$y]->[($totalgenes*2)+2];
			$all_diseases =~ s/;(0;)+/;/g;
			$all_diseases =~ s/;(1;)+/;/g;
			$all_diseases =~ s/\A0;//g;
			$all_diseases =~ s/\A1;//g;
			$all_diseases =~ s/;0\Z//g;
			$all_diseases =~ s/;1\Z//g;
			$all_diseases =~ s/\A1\Z/N\/A/g;
			$all_diseases =~ s/\A0\Z/N\/A/g;
			print OUT "$escribir"."\t".$all_diseases."\n";
		}
		close ARCHIVO3;
		close OUT;
		if ($s==1) {
			# unlink("$archivo1.input") or die "Could not delete the file!\n"; # To GitHub
			unlink("./$archivo1.input") or die "Could not delete the file!\n";
		}
		$loop++;
		if(($loop % 1000) == 0) {
			print "$loop analyzed\n";
		}
	}
}

# Closing Statements

close ARCHIVO1;
close ARCHIVO2;
close ARCHIVOG1;
close ARCHIVOG2;
close ARCHIVOG3;
close ARCHIVOG4;
close ARCHIVOG5;
close ARCHIVOG6;
close ARCHIVOG7;
print "\nDone!\n";