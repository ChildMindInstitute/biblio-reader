<?php

mysql_connect('localhost', 'root', 'renew');
mysql_select_db('annotate');

ob_end_clean();
ob_start();
echo "Debug output from PubMed search:\n\n";
chdir("/var/www/html/annotate");
passthru("python utils/rsSearch.py /var/www/html/annotate/working/add_these_new_pmids.txt");

exec("python utils/pmidToBib.py /var/www/html/annotate/working/add_these_new_pmids.txt");

exec("python utils/bibToSQL.py /var/www/html/annotate/working/add_these_new_pmids.bib");

passthru("mysql -u root --password=renew annotate < /var/www/html/annotate/working/add_these_new_pmids.sql");
$firstPart = ob_get_contents();
ob_end_clean();

ob_start();
echo "Annotator URL: http://ec2-54-224-192-157.compute-1.amazonaws.com/annotate/\n\n";
passthru("php utils/createExpertAssignments.php");
$secondPart = ob_get_contents();
ob_end_clean();

if (strpos($secondPart, "are being assigned") !== false) {
	$subject = 'New RS publications for review';	
} else {
	$subject = "No new RS publications";	
}

$headers = 'From: annotator@cmi-annotator.com' . "\r\n" .
    'Reply-To: annotator@cmi-annotator.com' . "\r\n" .
    'X-Mailer: PHP/' . phpversion();


require_once '/var/www/html/annotate/swift/lib/swift_required.php';
$transport = Swift_SmtpTransport::newInstance('smtp.gmail.com', 465, 'ssl')
  ->setUsername("matthew.k.doherty@gmail.com")
  ->setPassword("ii3q00rii3q00r");
$mailer = Swift_Mailer::newInstance($transport);
$message = Swift_Message::newInstance($subject)
  ->setFrom(array('matthew.k.doherty@gmail.com' => "Matt Doherty"))
  ->setTo(array('mattd@alum.mit.edu'))
  ->setBody($secondPart . "\n\n" . $firstPart)
  ;

// Send the message
$result = $mailer->send($message);  

?>