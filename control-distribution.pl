for(my $sgs=1; $sgs<=10; $sgs++){
my @ntgenes;
foreach my $i (1..5000000){
	my $score=0;
	for(my $x=0; $x<$sgs; $x++){
		$score+=log(rand(1));
	}
	$score/=$sgs;
	push @ntgenes, int($score*100)/100;
}
my @ntsort = sort {$a<=>$b} @ntgenes;
open(FILE, ">ctrlscores$sgs");
my $avg=0;
foreach my $score (@ntsort){
	$avg+=$score;
	$count++; 
	if($count==5){$avg/=$count; $count=0;print FILE "$avg\n";  $avg=0;}
}
if($count>0){$avg/=$count; $count=0;print FILE "$avg\n";  $avg=0;}
close FILE;
}
