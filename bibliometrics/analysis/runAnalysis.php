<?php
ob_implicit_flush();
ob_end_flush();
echo "<pre>\n";
chdir("/var/www/html/");
echo "Running overall growth...\n";
passthru("python growth/overallGrowth.py > figuresGen/working/overallGrowth.txt");
chdir("/var/www/html/figuresGen/");
echo "Plotting...\n";
passthru("Rscript biblio_helper.R");
echo "Done.\n\n<a href='/figures/'>View figures</a>\n";
echo "</pre>\n";
?>
