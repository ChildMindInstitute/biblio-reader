<?php
header('Content-type: text/html; charset=iso-8859-1');

mysql_connect('localhost', 'root', 'renew');
mysql_select_db('annotate');

$rs = mysql_query("select * from tags");
$tags = array();
while ($row = mysql_fetch_assoc($rs)) {
	$tags[] = $row;
}
foreach ($tags as &$tag) {
	$tagName = " " . mysql_real_escape_string(strtolower($tag['tags_name']));
	$rs = mysql_query("select refs.ref_id, assignment_id from refs, expert_assignments where refs.ref_id=expert_assignments.ref_id and (lower(abstract) like '%{$tagName}%' or lower(title) like '%{$tagName}%' or lower(mendeley_tags) like '%{$tagName}%')");
	$inserts = array();
	while ($row = mysql_fetch_assoc($rs)) {
		$inserts[] = "(" . $tag['tags_id'] . ", " . $row['assignment_id'] . ")";		
	}
	mysql_query("insert into tags_assignments (ta_tags_id, ta_assignments_id) values " . implode(",", $inserts));
}
