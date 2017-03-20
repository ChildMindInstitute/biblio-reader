<?php
header('Content-type: text/html; charset=iso-8859-1');

$tags = array('clinical','basic_neuroscience','brain_and_behavior','functional_anatomy','multimodal','genetics','animal_models','methodology','review_meta_analysis','1000_functional_connectomes','graph_theory',
			'supervised_learning','unsupervised_learning','independent_component_analysis','seed_based_correlation');
$internalTags = array('expert_needed', 'not_resting_state');
$clinicalTags = array();
$ignoreFields = array('ref_id', 'tags', 'assignment_id', 'master_tags', 'expert_name', 'assigned', 'author', 'affiliation', 'age_range_low', 'age_range_high', 'update_time', 'ta_id', 'ta_assignments_id', 'ta_tags_id');			
$maxResults = 50;

mysql_connect('localhost', 'root', 'renew');
mysql_select_db('annotate');

function getChosenTagsString($assId) {
	$ary = mysql_fetch_array(mysql_query("select group_concat(tags_name separator ',') from tags, tags_assignments where tags_id=ta_tags_id and ta_assignments_id=" . $assId));
	return $ary[0];
}

function prettyTag($tag) {
	return str_replace(' ', '&nbsp;', ucwords(str_replace('_', ' ', $tag)));
}

$rs = mysql_query("select tags_id, tags_name from tags where tags_hidden=0");
while ($ary = mysql_fetch_array($rs)) {
	$clinicalTags[$ary[0]] = $ary[1];
}

$tagDefs = array();
$tagDefs["Clinical"] = "Related to treatment of patients.";
$tagDefs["Basic Neuroscience"] = "Relevant to our understanding of the brain and nervous system.";
$tagDefs["Brain and Behavior"] = "Relates resting state to the behavior of subjects.";
$tagDefs["Functional Anatomy"] = "Uses resting state to describe the function of parts of the brain.";
$tagDefs["Genetics"] = "Relates resting state to genes and heredity.";
$tagDefs["Animal Models"] = "Uses animals to better understand humans.";
$tagDefs["Methodology"] = "Describes a novel method, such as a data processing algorithm.";
$tagDefs["Review/Meta-analysis"] = "Summarizes recent literature.";
$tagDefs["1000 Functional Connectomes"] = "Uses data provided by 1000 FC, a data-sharing initiative.";
$tagDefs["Multimodal"] = "Uses multiple imaging methods, such as fMRI and EEG.";
$tagDefs["Graph Theory"] = "Models data as a graph, with nodes and edges connecting them. For example, a node can represent an ROI and an edge its correlation with another ROI. A popular algorithm that operates on graphs is eigenvector centrality.";
$tagDefs["Supervised Learning"] = "Uses algorithms that build models from a labeled training set to make predictions about a test set. An example is the support vector machine (SVM).";
$tagDefs["Unsupervised Learning"] = "Uses algorithms that build models based on the structure of the data, without a labeled training set. An example is k-nearest neighbors (k-NN).";
$tagDefs["Independent Component Analysis"] = "A particular unsupervised learning algorithm that learns a small number of informative data vectors from a larger data set. Finds these vectors by maximizing their statistical independence. We have a separate tag for ICA because it's so useful in resting state analysis. Don't automatically tag ICA papers with the unsupervised learning tag.";
$tagDefs["Seed Based Correlation"] = "Measures connectivity between a seed voxel (or ROI) and the rest of the brain using correlation.";
$internalDefs = array();
$internalDefs["Expert Needed"] = "Set this to indicate you need help.";
$internalDefs["Not Resting State"] = "Our PubMed search is not as specific as it is sensitive. Use this for papers that aren't actually resting state.";
$internalDefs["Pending Approval"] = "Papers get this status after they're returned from the PubMed search. This means we're not yet sure whether this paper is about resting state.";
$internalDefs["Pending Tagging"] = "The paper is resting state, but it hasn't been completely tagged yet.";
$internalDefs["Tagging Complete"] = "The paper is resting state and tagging is complete. The expert's job is done.";
$internalDefs["Trash"] = "The paper is not resting state. The expert's job is done.";
$qa = array();



?>
<html>
	<head>
	<link rel="stylesheet" href="../figures/style.css" type="text/css" />
	<link rel="stylesheet" href="../figures/interior.css" type="text/css" />

	<style type="text/css">
	div
	{
		background-color: white;
	}
	.checkboxes label
	{
	    display: block;
	    float: center;
	    padding-right: 10px;
	    padding-bottom: 5px;
	    white-space: nowrap;
	}
	
	.checkboxes input
	{
	    vertical-align: middle;
	}
	
	.checkboxes label span
	{
	    vertical-align: middle;
	}
	div.content
	{
		text-align:left; width: 40%; display:inline-block;
	}
	p
	{
		padding-bottom: 20px;
	}
	</style>
	</head>    
	<body style="background-color: white;"> <!--"font:90%/200% Arial,Helvetica,sans-serif;"-->
		<div style="background-color: white; text-align:center; list-style-position:inside; font-size: 18px;">
			<p><h1>CMI Annotator Help</h1></p>
			<p>part of the <a href="http://www.childmind.org/en/healthy-brain-network/cmi-librarian/">CMI Librarian</a></p>
			<p>back to <a href="/annotate/index.php">Annotator</a></p>
			<p>&nbsp;</p>
			<p><b>Overview</b></p>
			<div class="content">
				<p>
					We currently (as of October) have several thousand resting state papers that need to be correctly tagged by our human experts. 
				</p>
				<p>
					The papers have been divided up into several sets, one for each expert. 20% of the papers are shared among all experts so we can get reliability statistics.					 
				</p>
				<p>
					To get started, click on the <a href="/annotate/">Annotator</a> link, then click your name. You'll then have access only to the papers that are assigned to you. Navigate through the papers using the options on the left, and tag them using the checkboxes. When you're finished with a paper, set its status to 'Tagging Complete' or 'Trash'. Statuses for all the experts are summarized on the <a href="/figures/">figures</a> page.
				</p>
				<p>
					The papers' tags were initialized with guesses that I made using simple text matching just to get the ball rolling, so don't assume that tags that are already there under your name are correct.
				</p>
			</div>
			<p>&nbsp;</p>
			<p><b>Tag Definitions</b></p>
			<div class="content">
				<?php
				foreach ($tagDefs as $key => $value) {
					echo "<p><b>{$key}</b> - {$value}</p>";
				}
				?>
			</div>
			<p>&nbsp;</p>
			<p><b>Internal Definitions</b></p>
			<div class="content">
				<?php
				foreach ($internalDefs as $key => $value) {
					echo "<p><b>{$key}</b> - {$value}</p>";
				}
				?>
			</div>
			<p>&nbsp;</p>
			<p><b>Q&A</b></p>
			<div class="content">
				<?php
				foreach ($qa as $key => $value) {
					echo "<p><b>{$key}</b> - {$value}</p>";
				}
				?>
			</div>				
		</div>
	</body>
</html>
