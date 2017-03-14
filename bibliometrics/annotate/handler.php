<?php

print_r($_REQUEST);
if (isset($_REQUEST['ref_id'])) {
	$refId = $_REQUEST['ref_id'];
	$tags = "";
	$ageRangeLow = "";
	$ageRangeHigh = "";
	$sampleSize = "";
	mysql_connect('localhost', 'root', 'renew');
	mysql_select_db('annotate');
	/*
	if (isset($_REQUEST['tags'])) {
		$tags = implode(",", $_REQUEST['tags']);
	}*/
	if (isset($_REQUEST['tags'])) {
		$tags = $_REQUEST['tags'];
	}
	if (isset($_REQUEST['ageRangeLow'])) {
		$ageRangeLow = $_REQUEST['ageRangeLow'];
	}
	if (isset($_REQUEST['ageRangeHigh'])) {
		$ageRangeHigh = $_REQUEST['ageRangeHigh'];
	}
	if (isset($_REQUEST['sampleSize'])) {
		$sampleSize = $_REQUEST['sampleSize'];
	}
	/*
	if (isset($_REQUEST['expert_name'])) {
		$expertName = $_REQUEST['expert_name'];
		$query = 'update expert_assignments set tags="' . $tags . '" where ref_id="' . $refId . '" and expert_name="' . $expertName . '" limit 1';
		echo $query;
		$result = mysql_query($query);
		
	} else {
		die("NO EXPERT");
	}	
	 */
	if (isset($_REQUEST['assignment_id'])) {
		$assId = $_REQUEST['assignment_id'];
		mysql_query("delete from tags_assignments where ta_assignments_id=" . $assId);
		$checkedParents = array();
		foreach ($tags as &$tagId) {
			mysql_query("insert into tags_assignments (ta_tags_id, ta_assignments_id) values (" . $tagId . "," . $assId . ")");
			echo "inserted " . $tagId . "\n";
			$rs = mysql_query("select tags_parent_id from tags where tags_id={$tagId}");
			$row = mysql_fetch_assoc($rs);
			$parentId = $row['tags_parent_id'];
			if ($parentId > 0 && in_array($parentId, $checkedParents) === false && in_array($parentId, $tags) === false) {
				$checkedParents[] = $parentId;
			}
		}
		$query = "update expert_assignments set age_range_low='{$ageRangeLow}', age_range_high='{$ageRangeHigh}', sample_size='{$sampleSize}' where assignment_id={$assId}";
		echo $query;
		mysql_query($query);
		echo "\nCHECKED_PARENTS:" . implode(",", $checkedParents);
	} else {
		die("NO EXPERT");
	}	
	
}
?>
