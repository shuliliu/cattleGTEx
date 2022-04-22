$geno_vcf = $ARGV[0];

$query_bed = $ARGV[1];





open FILE, "$ARGV[3]"; #covariance not gzipped

while(<FILE>){

   next if $_ !~ /\d/;

     

   my @data = split /\s+/, $_;

   shift @data until $data[0]=~/^\S/;



   if (/\s*id/){

      shift @data;

      @id_list = @data;

      next;

   }



   $var = $data[0];

   shift @data;

   $covar_text .="controlled $var gtex @data\n";

}

















$geno_text = "";

open FILE, "bcftools view -R $query_bed $geno_vcf |"; #mean the vcf file should be gzipped and tabix

while(<FILE>){



   if(/^\s*\#CHROM/){

      my @data = split /\s+/,$_;



      shift @data until $data[0]=~/^\S/;
      shift @data for (1..9);

      @hash{@data} = (0..$#data);

      foreach $id (@id_list){

         push @index, $hash{$id};

      }

      next;

   }

   $geno_text.="@id_list\n";
}
   print "$geno_text";