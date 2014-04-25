<?php
header('Content-type: text/html; charset=iso-8859-1');

mysql_connect('localhost', 'root', 'renew');
mysql_select_db('annotate');
$overlapFraction = 0.2;
$numExperts = 7;

$overallCountAry = mysql_fetch_array(mysql_query('select count(*) from refs where assigned=0'));
$overallCount = $overallCountAry[0];
echo "Number of new publications: " . $overallCount . "\n";

if ($overallCount == 0) {
	echo "There are no new resting state publications for review\n\n";
	exit();
}

echo "The following publications are being assigned to experts for review with status 'pending approval':\n\n";
$result = mysql_query('select title, author, pmid from refs where assigned=0 order by title');
while ($ary = mysql_fetch_assoc($result)) {
	echo "Title: " . mb_convert_encoding($ary['title'], "ASCII", "UTF-8") . "\n";
	echo "Author: " . mb_convert_encoding($ary['author'], "ASCII", "UTF-8") . "\n";
	echo "PMID: " . $ary['pmid'] . "\n";
	echo "\n";
}

echo "Debug output from assigner:\n";
$overlapCount = intval(($overallCount / $numExperts) * $overlapFraction);
echo "\nOverlap count: " . $overlapCount;

$start = $overlapCount;
$countEach = intval(($overallCount - $overlapCount) / $numExperts + 1);
echo "\nCount each: " . $countEach;

for ($i = 0; $i < $numExperts; $i++) {
	$result = mysql_query('select ref_id from refs where assigned=0 order by title limit ' . $start . ', ' . $countEach);
	// XXX???
	//$result = mysql_query('select refs.ref_id as ref_id, expert_assignments.tags as tags, expert_assignments.status as status from refs, expert_assignments where refs.ref_id=expert_assignments.ref_id and assigned=0 and expert_name="Master" order by title limit ' . $start . ', ' . $countEach);
	while ($ary = mysql_fetch_assoc($result)) {
		$tuple = implode(",", array("'Expert" . ($i+1) . "'", $ary['ref_id'], "'" . $ary['tags'] . "'", "'" . $ary['status'] . "'"));
		mysql_query('insert into expert_assignments (expert_name, ref_id, tags, status) values (' . $tuple . ')');
	}
	echo "\nFor Expert" . $i . " we inserted: " . $start . " to " . ($start+$countEach-1);
	$start += $countEach;	
}

//$result = mysql_query('select ref_id from refs where assigned=0 order by title limit 0, ' . $overlapCount);
$result = mysql_query('select refs.ref_id, expert_assignments.tags, expert_assignments.status from refs, expert_assignments where refs.ref_id=expert_assignments.ref_id and assigned=0 and expert_name="Master" order by title limit 0, ' . $overlapCount);
while ($ary = mysql_fetch_assoc($result)) {
	for ($i = 0; $i < $numExperts; $i++) {
		$tuple = implode(",", array("'Expert" . ($i+1) . "'", $ary['ref_id'], "'" . $ary['tags'] . "'", "'" . $ary['status'] . "'"));
		mysql_query('insert into expert_assignments (expert_name, ref_id, tags, status) values (' . $tuple . ')');
	}
}
echo "\nFor all experts inserted: 0 to " . $overlapCount;

$result = mysql_query('select ref_id, master_tags from refs where assigned=0 order by title');
while ($ary = mysql_fetch_assoc($result)) {
	//mysql_query('insert into expert_assignments (expert_name, ref_id, tags) values ("Master", ' . $ary['ref_id'] . ', "' . $ary['master_tags'] . '")');
}
//echo "\nFor master inserted " . $overallCount . "\n";

mysql_query("update refs set assigned=1 where assigned=0");
?>