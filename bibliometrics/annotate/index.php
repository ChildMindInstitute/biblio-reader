<?php
ini_set('display_errors',1);
ini_set('display_startup_errors',1);
error_reporting(-1);
header('Content-type: text/html; charset=iso-8859-1');

$tags = array('clinical','basic_neuroscience','brain_and_behavior','functional_anatomy','multimodal','genetics','animal_models','methodology','review_meta_analysis','1000_functional_connectomes','graph_theory',
			'supervised_learning','unsupervised_learning','independent_component_analysis','seed_based_correlation');
$internalTags = array('expert_needed', 'not_resting_state');
$tags = array();
$clinicalTags = array();
$ignoreFields = array('ref_id', 'tags', 'assignment_id', 'master_tags', 'expert_name', 'assigned', 'author', 'affiliation', 'age_range_low', 'age_range_high', 'sample_size', 'update_time', 'ta_id', 'ta_assignments_id', 'ta_tags_id', 'volume', 'pages', 'number', 'isbn', 'journal_iso');			
$maxResults = 10;

mysql_connect('localhost', 'root', 'renew');
mysql_select_db('annotate');

function getChosenTagsString($assId) {
	$ary = mysql_fetch_array(mysql_query("select group_concat(tags_name separator ',') from tags, tags_assignments where tags_id=ta_tags_id and ta_assignments_id=" . $assId));
	return $ary[0];
}

function prettyTag($tag) {
	if ($tag == "review_meta_analysis") {
		return "Review/Meta-analysis";
	} else if ($tag == "doi") {
		return "DOI";
	} else if ($tag == "pmid") {
		return "PMID";		
	} else {
		return str_replace(' ', '&nbsp;', ucwords(str_replace('_', ' ', $tag)));
	}
}

function echoTagHtml($tagId, $tags, $chosenTags, $level=0) {
	$tagName = $tags[$tagId];
?>
<label style='margin-left:<?php echo ($level*10) . "px;" ?><?php if (in_array($tagName, explode(",", $chosenTags)) !== false) { echo "font-weight: bold;"; } ?>'><input id="tag<?php echo $tagId; ?>" type='checkbox' name='tags[]' <?php if (in_array($tagName, explode(",", $chosenTags)) !== false) { echo "checked"; } ?> value='<?php echo $tagId; ?>'>&nbsp;<span><?php echo $tagName; ?></span></label>
<?php		
}

function getCheckForm($nextIx) {
	return "<form id='checks' action='index.php#aid" . $nextIx . "' method='post'>";
}

function getExpertLinks() {
	global $expert_name;
?>			
			<p>
				<?php 
				$result = mysql_query("select distinct expert_name from expert_assignments");
				while ($ary = mysql_fetch_array($result)) {					
					if ($expert_name == $ary['expert_name']) {
						echo $expert_name . "&nbsp; ";
					} else {
						echo "<a href='index.php?expert_name=" . $ary['expert_name'] . "'>" . $ary['expert_name'] . "</a>&nbsp; "; 
					}
				}
				echo "<br><a href='/figures'>View Figures</a>&nbsp; <a href='help.php'>Help</a>&nbsp; <a href='https://cmi.hackpad.com/CMI-Annotator-PRbcvRqgaxj'>HackPad Notes</a></b>";
				?>
			</p>
<?php
}

$expert_name = "";
if (isset($_REQUEST['expert_name'])) {
	$expert_name = $_REQUEST['expert_name'];
}

if (isset($_REQUEST['addClinicalTagSubmit']) && $_REQUEST['addClinicalTagSubmit']) {
	$newTag = mysql_real_escape_string($_REQUEST['addClinicalTag']);
	$rs = mysql_query("select tags_id, tags_hidden from tags where tags_name='{$newTag}' and tags_parent_id='20' limit 1");
	if (mysql_num_rows($rs) > 0) {
		$row = mysql_fetch_assoc($rs);
		if ($row['tags_hidden'] == 1) {
			mysql_query("update tags set tags_hidden=0 where tags_id={$row['tags_id']} and tags_parent_id='20' limit 1");
		}
	} else {
		mysql_query("insert into tags (tags_name, tags_parent_id) values ('{$newTag}', '20')");
	}
}

if (isset($_REQUEST['deleteClinicalTagSubmit']) && $_REQUEST['deleteClinicalTagSubmit']) {
	$deleteTagId = $_REQUEST['clinicalTag'];
	mysql_query("update tags set tags_hidden=1 where tags_id={$deleteTagId} and tags_parent_id='20' limit 1");
}

$rs = mysql_query("select tags_id, tags_name from tags where tags_hidden=0 and tags_parent_id=20 order by tags_name");
while ($ary = mysql_fetch_array($rs)) {
	$clinicalTags[$ary[0]] = $ary[1];
}

$rs = mysql_query("select tags_id, tags_name from tags where tags_hidden=0 and tags_parent_id!=20 order by tags_name");
while ($ary = mysql_fetch_array($rs)) {
	$tags[$ary[0]] = $ary[1];
}


if (isset($_REQUEST['mark_refresh_button']) && $_REQUEST['mark_refresh_button']) {
	$val = $_REQUEST['mark_refresh_button'];
	if (isset($_REQUEST['ref_id'])) {
		$refId = $_REQUEST['ref_id'];
		$tagAry = array();
		if (isset($_REQUEST['tags'])) {
			$tagAry = $_REQUEST['tags'];
		}
		if ($val == "Mark PENDING TAGGING") {			
			$newStatus = "pending_tagging";
		} else if ($val == "Mark APPROVED AND PENDING TAGGING") {			
			$newStatus = "pending_tagging";
		} else if ($val == "Mark REMOVE THIS") {			
			$newStatus = "trash";
		} else if ($val == "Mark PENDING APPROVAL") {			
			$newStatus = "pending_approval";
		} else if ($val == "Mark TAGGING COMPLETE") {			
			$newStatus = "tagging_complete";
		} else {
			die("Val: " . $val);
		}

		mysql_connect('localhost', 'root', 'renew');
		mysql_select_db('annotate');
		if (isset($_REQUEST['expert_name'])) {
			$expertName = $_REQUEST['expert_name'];
			$query = 'update expert_assignments set status="' . $newStatus . '" where ref_id="' . $refId . '" and expert_name="' . $expertName . '" limit 1';
		} else {
			die("NO EXPERT");
			$query = 'update refs set master_tags="' . $tagUpdate . '" where ref_id="' . $refId . '" limit 1';
		}	
		$result = mysql_query($query);
	}
}

$completeTags = array();
$break = 0;
if (isset($_REQUEST['break'])) {
	$break = $_REQUEST['break'];
}

$extraWhere = "";
$extraTables = "";
$addedTa = false;
if (isset($_REQUEST['tag']) && $_REQUEST['tag'] != "") {
	$tag = $_REQUEST['tag'];
	$completeTags[] = prettyTag($tags[$tag]);
	$extraWhere .= " and ta_tags_id={$tag} and ta_assignments_id=assignment_id";
	$extraTables .= ', tags_assignments';
	$addedTa = true;
}
if (isset($_REQUEST['clinicalTag']) && $_REQUEST['clinicalTag'] != "") {
	$clinicalTag = $_REQUEST['clinicalTag'];
	$completeTags[] = prettyTag($clinicalTags[$clinicalTag]);
	$extraWhere .= " and ta_tags_id={$clinicalTag} and ta_assignments_id=assignment_id";
	if (!$addedTa) {
		$extraTables .= ', tags_assignments';
	}
	$addedTa = true;
}
if (isset($_REQUEST['year']) && $_REQUEST['year'] != "") {	
	$year = $_REQUEST['year'];
	$extraWhere .= ' and year="' . $year . '" ';
}
if (isset($_REQUEST['freetext']) && $_REQUEST['freetext'] != "") {	
	$freetext = mysql_real_escape_string($_REQUEST['freetext']);
	$extraWhere .= " and (title like '%{$freetext}%' or abstract like '%{$freetext}%' or pmid='{$freetext}')";
}
if (isset($_REQUEST['rsFmri']) && $_REQUEST['rsFmri'] != "") {	
	$rsFmri = $_REQUEST['rsFmri'];
	$extraWhere .= " and (not exists (select * from tags_assignments where (ta_tags_id=55 or ta_tags_id=56) and ta_assignments_id=assignment_id))";
	if ($addedTa) {
		
	} else {
		
	}
}

$status = "pending_tagging";
if (isset($_REQUEST['status']) && $_REQUEST['status'] != "") {
	$status = $_REQUEST['status'];
}
if ($status != "any") {
	$extraWhere .= ' and status="' . $status . '" ';
} 

$refs = array();
$breaks = array();
$overallCount = 0;
if ($expert_name) {
	$mainQuery = 'from refs, expert_assignments' . $extraTables . ' where refs.ref_id = expert_assignments.ref_id and expert_assignments.expert_name="' . $expert_name . '" ' . $extraWhere . ' order by title';
	$result = mysql_query("select refs.ref_id as ref_id, title, author, pmid, doi, journal, year, volume, pages, number, abstract, affiliation, isbn, journal_iso, mesh_heading_list, assigned, mendeley_tags, source, assignment_id, tags, status, age_range_low, age_range_high, sample_size, update_time {$mainQuery} limit {$break}, {$maxResults}") or die(mysql_error());
	$overallCountAry = mysql_fetch_array(mysql_query("select count(*) {$mainQuery}"));
	$overallCount = $overallCountAry[0];
	
	$completeTagText = implode(', ', $completeTags);
	
	while ($row = mysql_fetch_assoc($result)) {
		$refs[] = $row;
	}
	
	for ($i = 0; $i < $overallCount; $i += $maxResults) {
		$breaks[] = $i;
	}
}

?>

<!DOCTYPE html>
<!-- saved from url=(0064)http://www.childmind.org/en/healthy-brain-network/cmi-librarian/ -->
<html lang="en" class="wf-locatorweb1-n7-active wf-locatorweb1-n2-active wf-locatorweb1-n4-active wf-active"><!--<![endif]--><head><meta http-equiv="Content-Type" content="text/html; charset=UTF-8"><script type="text/javascript" src="./cmi_librarian_files/shares.json"></script><script src="./cmi_librarian_files/cb=gapi.loaded_1" async=""></script><script src="./cmi_librarian_files/cb=gapi.loaded_0" async=""></script><script type="text/javascript" src="./cmi_librarian_files/auth016.js"></script><script id="LR1" type="text/javascript" async="" src="./cmi_librarian_files/client.js"></script><script type="text/javascript" src="./cmi_librarian_files/counter017.js"></script><link rel="stylesheet" type="text/css" href="./cmi_librarian_files/counter014.css" media="all"><script type="text/javascript" src="./cmi_librarian_files/plusone.js" gapi_processed="true"></script><script type="text/javascript" src="./cmi_librarian_files/widgets.js"></script><link rel="stylesheet" type="text/css" href="./cmi_librarian_files/widget119.css" media="all"><script src="./cmi_librarian_files/addthis_widget.js"></script><script src="./cmi_librarian_files/global.js"></script><script src="./cmi_librarian_files/flexie.js"></script><script src="./cmi_librarian_files/jquery.cookie.js"></script><script src="./cmi_librarian_files/jquery.easing.js"></script><script src="./cmi_librarian_files/jquery-ui-1.8.7.custom.min.js"></script><script src="./cmi_librarian_files/jquery-1.4.2.js"></script>
		<meta charset="utf-8">

		<!-- http://www.phpied.com/conditional-comments-block-downloads/ -->
		<!--[if IE]><![endif]-->

		<meta http-equiv="X-UA-Compatible" content="IE=edge,chrome=1">

		<title> - CMI Librarian  | Child Mind Institute</title>

		<meta name="description" content="The Child Mind Institute (CMI) Librarian initiative aims to facilitate the aggregation, distillation and dissemination of findings and insights, through the creation and sharing of hand-vetted and sorted reference libraries.">
		<meta name="author" content="Child Mind Institute">

		<meta name="viewport" content="width=device-width">

		<link rel="shortcut icon" href="http://dafjmlate1fc5.cloudfront.net/favicon.ico">
		<link rel="apple-touch-icon" href="http://dafjmlate1fc5.cloudfront.net/apple-touch-icon.png">

		<!--[if ! lte IE 6]><!-->
		<!-- Global CSS -->
		<link rel="stylesheet" href="./cmi_librarian_files/style.css" type="text/css">

		<!-- Local CSS -->
		<link rel="stylesheet" href="./cmi_librarian_files/interior.css" type="text/css">

		<!--<script id="twitter-wjs" src="./cmi_librarian_files/widgets.js"></script><script type="text/javascript" async="" src="./cmi_librarian_files/ga.js"></script><script type="text/javascript" src="./cmi_librarian_files/usq2wty.js"></script><style type="text/css"></style>-->
		<!--<style type="text/css">.tk-locator-web{font-family:"locator-web-1","locator-web-2",sans-serif;}</style><link rel="stylesheet" href="http://use.typekit.com/c/8ff65e/locator-web-1:n2:n4:n7.PZ1:M:2,PYx:M:2,PYz:M:2/d?3bb2a6e53c9684ffdc9a9bff1f5b2a626519ae8d4e1a1ed3b663be3f1c043f7b52411b0aca68ee6ca721baf897a2004d81a9e919cc964b46e1093a7abe96fbf23ae8da5760bac5b5d4725282644756f1df338554bab51510cd37804f06ea787c9775c232d4121b334a3f9b732c6d10fb5de96517e04c692e4effdf043f9218fcf47bd35516d06a09fff4885c08c83c75eb4f26dff2a29beb6893f857682aca5a1f672b6e2b72f5fe1727e369ce429d3a7f4f8ec464d4572f6f411f17a869562c8af02c"><script type="text/javascript">try{Typekit.load();}catch(e){}</script>-->
		<!--<![endif]-->

		<!--[if lte IE 6]>
		<link rel="stylesheet" href="/static/css/ie6.css?v=1.8.4" media="screen, projection" />
		<![endif]-->

		<!--[if lt IE 9]>
		<script src="/static/js/cmn/lib/plugins/IE9.js?v=1.8.4"></script>
		<![endif]-->
	<style type="text/css">

	body #background-shim {
		background: #dbdcdd url("cmi_librarian_files/bg-stripe.png") repeat-x;	
		position: absolute;	
	}

	html {
		background: #ecedef url("cmi_librarian_files/bg.png") repeat-y top;
		height: 100%;
	}

	.column.primary .heading {
		background: #024990 url("cmi_librarian_files/bg-heading.jpg") no-repeat;
	}
	footer {	
		background: #ecedef url("cmi_librarian_files/bg-footer.png") no-repeat center top;
	}
	div
	{
		#background-color: white;
	}
	.checkboxes label
	{
	    display: block;
	    float: center;
	    padding-right: 10px;
	    padding-top: 5px;
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

	form
	{
		margin: 0px;
		padding: 0px;
	}
	.anno_heading {
		padding-top: 20px;
	}
	td {
		padding-top: 5px;
		padding-bottom: 5px;
		padding-left: 5px;
		line-height: 150%;
	}
	
	input[type='text'],section.module form input[type='password'],section.module form select{background-color: white; clear:left;padding: 4px;border:1px solid #dbdbdb;-webkit-box-shadow:rgba(0,0,0,0.15),0,3px,3px,0,inset;-moz-box-shadow:rgba(0,0,0,0.15),0,3px,3px,0,inset;box-shadow:rgba(0,0,0,0.15),0,3px,3px,0,inset}	
	</style>
	<script type="text/javascript" src="http://code.jquery.com/jquery-1.10.0.js"></script>    
    <script type="text/javascript" charset="utf-8">
      // let's start the jQuery while I wait.
      // step 1: onload - capture the submit event on the form.
      $(function() { // onload...do
      	var request;
        var doUpdate = function() {
		    // abort any pending request
		    if (request) {
		        request.abort();
		    }
		    // setup some local variables
		    var $form = $(this).closest('form');
		    // let's select and cache all the fields
		    var $inputs = $form.find("input, select, button, textarea");
		    // serialize the data in the form
		    var serializedData = $form.serialize();
			console.log("Serialized data: " + serializedData);
		    // let's disable the inputs for the duration of the ajax request
		    //$inputs.prop("disabled", true);
		
		    // fire off the request to /form.php
		    request = $.ajax({
		        url: "handler.php",
		        type: "post",
		        data: serializedData
		    });
		
		    // callback handler that will be called on success
		    request.done(function (response, textStatus, jqXHR){
		        // log a message to the console
		        console.log("Hooray, it worked!" + response);
		        var marker = "CHECKED_PARENTS:";
		        var subby = response.substr(response.indexOf(marker)+marker.length);
		        var checkedParents = subby.split(",");
		        checkedParents.forEach(function(parent) {       			        	
		        	if (parent > 0) {
		        		var thing = $form.find('#tag' + parent)
			        	console.log("Sending click to " + parent + " " + thing.attr('name'));
			        	thing.click();
			        }        	
		        });
		    });
		
		    // callback handler that will be called on failure
		    request.fail(function (jqXHR, textStatus, errorThrown){
		        // log the error to the console
		        console.error(
		            "The following error occured: "+
		            textStatus, errorThrown
		        );
		    });
		
		    // callback handler that will be called regardless
		    // if the request failed or succeeded
		    request.always(function () {
		        // reenable the inputs
		        $inputs.prop("disabled", false);
		    });
		
		    // prevent default posting of form
		    //event.preventDefault();
            // re-test...
            // by default - we'll always return false so it doesn't redirect the user.          
		    if ($(this).is(':checked')) $(this).closest('label').css('font-weight', 'bold');
            else $(this).closest('label').css('font-weight', 'normal');
            return true;
		}
		$(document).on('click', ':checkbox', doUpdate);
		$(document).on('input', '[id=ageLow]', doUpdate);
		$(document).on('input', '[id=ageHigh]', doUpdate);
		$(document).on('input', '[id=sampleSize]', doUpdate);
      })
    </script>

</head>

	<body class="page">
		<!--[if ! lte IE 6]><!-->
		<script src="./cmi_librarian_files/bootstrap.js" type="text/javascript"></script>

		<!--<![endif]-->
		<div id="wrapper">
			
			<header>
				<div id="logo">
					<a href="http://www.childmind.org/"><img src="./cmi_librarian_files/logo.png" width="274" height="83" alt="Child Mind Institute Logo"></a>
				</div>

				<div id="superfluous-info">

					<div id="account-info">
						<ul>
							<li><a href="http://www.childmind.org/request-an-appointment/">Request an Appointment</a></li>
							<li><a href="http://www.childmind.org/en/announcing-patient-portal/">Patient Portal</a></li>
							<li>
								<ul class="social-icons">
									<li><a href="http://www.facebook.com/ChildMindInstitute" target="_blank"><img src="./cmi_librarian_files/facebook.png" alt="Facebook" width="20" height="20"></a></li>
									<li><a href="http://twitter.com/ChildMindDotOrg" target="_blank"><img src="./cmi_librarian_files/twitter.png" alt="Twitter" width="20" height="20"></a></li>
									<li><a href="https://plus.google.com/u/0/b/104916968679971130579/" target="_blank"><img src="./cmi_librarian_files/gplus.png" alt="Google+" width="20" height="20"></a></li>
								</ul>
							</li>
							<li>
							
							<a class="actionable-button modal" href="http://www.childmind.org/en/healthy-brain-network/cmi-librarian/#modal_login-register">Sign In</a>
							
							</li>
						</ul>
					</div>

					<div id="donation-panel">
						<ul>
							<li>
								<form action="http://www.childmind.org/search/" method="get">
									<fieldset>
										<input type="text" name="q" id="q" placeholder="Search" value="">
										<input type="submit" name="s" id="s" value="Search">
									</fieldset>
								</form>
							</li>
							<li>
								<a id="email-sign-up" href="http://www.childmind.org/en/healthy-brain-network/cmi-librarian/#email-sign-up" class="actionable-button large">Email Sign Up</a>
							</li>
							<li>
								<a href="https://secure3.convio.net/cmi/site/Donation2?1820.donation=form1&df_id=1820" class="actionable-button large callout" target="_blank">Donate</a>
							</li>
						</ul>
						<div id="email-overlay">
							<form method="post" action="http://support.childmind.org/site/Survey">
								<fieldset>
									<input type="hidden" name="SURVEY_ID" id="SURVEY_ID" value="1162">
									<input type="hidden" name="cons_info_component" id="cons_info_component" value="t">
									<input type="hidden" name="3726_1162_2_1622" id="3726_1162_2_1622_1" value="1044">
									<input type="hidden" name="3726_1162_2_1622" id="3726_1162_2_1622_2" value="1043">
									<input type="hidden" name="3726_1162_2_1622" id="3726_1162_2_1622_3" value="1361">
									<input type="text" name="cons_email" id="email" placeholder="Email Address">
									<input type="text" name="cons_zip_code" id="postal" size="10" placeholder="Zip Code">
									<input type="submit" name="ACTION_SUBMIT_SURVEY_RESPONSE" value="Submit">
									<input type="reset" value="Cancel">
								</fieldset>
							</form>
						</div>
					</div>

				</div>

				<nav>
					



<ul>
	
	
	<li class="get-information ">
		<a href="http://www.childmind.org/get-information/">Get Information</a>
        
        
		<ul>
		
			
			<li class="mental-health-guide ">
			<a href="http://www.childmind.org/health/disorder-guide/">Mental Health Guide</a></li>
			
		
			
			<li class="symptom-checker ">
			<a href="http://www.childmind.org/health/symptom-checker">Symptom Checker</a></li>
			
		
			
		
			
			<li class="developmental-milestones ">
			<a href="http://www.childmind.org/developmental-milestones/">Developmental Milestones</a></li>
			
		
			
			<li class="quick-facts ">
			<a href="http://www.childmind.org/quick-facts">Quick Facts</a></li>
			
		
			
			<li class="glossary ">
			<a href="http://www.childmind.org/glossary/">Glossary</a></li>
			
		
			
			<li class="other-resources ">
			<a href="http://www.childmind.org/other-resources/">Other Resources</a></li>
			
		
			
			<li class="parents-guide-to-getting-good-care ">
			<a href="http://www.childmind.org/parents-guide-getting-good-care/">Parents Guide to Getting Good Care</a></li>
			
		
			
		
			
		
		</ul>
		
        
	</li>
	
	
	
	<li class="advice-support ">
		<a href="http://www.childmind.org/advice-support/">Advice &amp; Support</a>
        
        
		<ul>
		
			
		
			
		
			
			<li class="hot-topics ">
			<a href="http://www.childmind.org/hot-topics/">Hot Topics</a></li>
			
		
			
			<li class="ask-a-clinician ">
			<a href="http://www.childmind.org/posts/ask-an-expert/">Ask an Expert</a></li>
			
		
			
			<li class="real-stories ">
			<a href="http://www.childmind.org/posts/real-stories/">Real Stories</a></li>
			
		
			
			<li class="in-the-news ">
			<a href="http://www.childmind.org/press/brainstorm/">Brainstorm Blog</a></li>
			
		
			
			<li class="free-workshops ">
			<a href="http://www.childmind.org/workshop-series/">Free Workshops</a></li>
			
		
			
		
			
		
			
		
		</ul>
		
        
	</li>
	
	
	
	<li class="find-treatment ">
		<a href="http://www.childmind.org/find-treatment/">Find Treatment</a>
        
        
		<ul>
		
			
		
			
			<li class="our-care ">
			<a href="http://www.childmind.org/our-care/">Our Care</a></li>
			
		
			
			<li class="our-approach ">
			<a href="http://www.childmind.org/our-approach/">Our Approach</a></li>
			
		
			
			<li class="clinician-staff-directory ">
			<a href="http://www.childmind.org/directory/clinicians/">Clinician Directory</a></li>
			
		
			
			<li class="what-to-expect ">
			<a href="http://www.childmind.org/what-to-expect/">What to Expect</a></li>
			
		
			
			<li class="request-an-appointment ">
			<a href="http://www.childmind.org/request-an-appointment/">Request an Appointment</a></li>
			
		
			
			<li class="who-else-can-help ">
			<a href="http://www.childmind.org/who-else-can-help/">Who Else Can Help?</a></li>
			
		
		</ul>
		
        
	</li>
	
	
	
	<li class="science-innovation current">
		<a href="http://www.childmind.org/science-innovation/">Science &amp; Innovation</a>
        
        
		<ul>
		
			
		
			
		
			
		
			
			<li class="center-for-the-developing-brain ">
			<a href="http://www.childmind.org/center-for-developing-brain/">Center for the Developing Brain</a></li>
			
		
			
		
			
			<li class="healthy-brain-network ">
			<a href="http://www.childmind.org/healthy-brain-network/">Healthy Brain Network</a></li>
			
		
			
		
			
		
			
			<li class="clinical-innovation ">
			<a href="http://www.childmind.org/clinical-innovation">Clinical Innovation</a></li>
			
		
			
			<li class="science-team-directory ">
			<a href="http://www.childmind.org/directory/scientists/">Science Team Directory</a></li>
			
		
			
			<li class="scientific-research-council ">
			<a href="http://www.childmind.org/scientific-research-council/">Scientific Research Council</a></li>
			
		
			
		
			
			<li class="distinguished-scientist-award ">
			<a href="http://www.childmind.org/distinguished-scientist-award-2014/">Distinguished Scientist Award</a></li>
			
		
			
		
			
			<li class="rising-scientist-award ">
			<a href="http://www.childmind.org/rising-scientist-award/">Rising Scientist Award</a></li>
			
		
		</ul>
		
        
	</li>
	
	
	
	<li class="get-involved ">
		<a href="http://www.childmind.org/get-involved/">Get Involved</a>
        
        
		<ul>
		
			
			<li class="give ">
			<a href="http://www.childmind.org/give/">Give</a></li>
			
		
			
			<li class="get-connected ">
			<a href="http://www.childmind.org/get-connected/">Get Connected</a></li>
			
		
			
			<li class="jobs-volunteering ">
			<a href="http://www.childmind.org/jobs-volunteering/">Jobs &amp; Volunteering</a></li>
			
		
			
			<li class="art-at-cmi ">
			<a href="http://www.childmind.org/art-at-cmi/">Student Art Project</a></li>
			
		
			
			<li class="working-at-the-child-mind-institute ">
			<a href="http://www.childmind.org/events/">Events</a></li>
			
		
			
		
			
			<li class="shop ">
			<a href="http://www.childmind.org/shop/">Shop</a></li>
			
		
			
			<li class="this-may-speak-up-for-kids ">
			<a href="http://speakup.childmind.org/">Speak Up for Kids</a></li>
			
		
		</ul>
		
        
	</li>
	
	
	
	<li class="about-us ">
		<a href="http://www.childmind.org/about-us/">About Us</a>
        
        
		<ul>
		
			
			<li class="leadership ">
			<a href="http://www.childmind.org/directory/leadership/">Executive Team</a></li>
			
		
			
			<li class="board-of-directors ">
			<a href="http://www.childmind.org/board-of-directors/">Board of Directors</a></li>
			
		
			
			<li class="our-scientific-research-council ">
			<a href="http://www.childmind.org/our-scientific-research-council/">Scientific Research Council</a></li>
			
		
			
			<li class="clinician-directory ">
			<a href="http://www.childmind.org/directory/clinicians/">Clinician Directory</a></li>
			
		
			
			<li class="scientists ">
			<a href="http://www.childmind.org/directory/scientists/">Science Team Directory</a></li>
			
		
			
			<li class="directory ">
			<a href="http://www.childmind.org/directory/">General Directory</a></li>
			
		
			
			<li class="test-tables ">
			<a href="http://www.childmind.org/our-partners/">Our Partners</a></li>
			
		
			
			<li class="annual-report ">
			<a href="http://www.childmind.org/annual-reports/">Annual Report</a></li>
			
		
			
			<li class="take-a-virtual-tour ">
			<a href="http://www.childmind.org/virtual-tour/">Take a Virtual Tour</a></li>
			
		
			
			<li class="press-room ">
			<a href="http://www.childmind.org/press/">Press Room</a></li>
			
		
			
			<li class="contact-us ">
			<a href="http://www.childmind.org/contact/">Contact Us</a></li>
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
			
		
		</ul>
		
        
	</li>
	
	
	
	
</ul>




				</nav>
			</header>

			<div id="content">

				

	<div class="column secondary">
		<h1><a href="http://www.childmind.org/en/healthy-brain-network/cmi-librarian/">CMI Librarian</a></h1>
		<nav>
        
		
<ul>
	
	<li class="healthy-brain-network current">
		<a href="index.php">Annotator</a>
        
        
		<ul>
		
				<?php 
				$result = mysql_query("select distinct expert_name from expert_assignments");
				while ($ary = mysql_fetch_array($result)) {
					?>
			<li class="the-child-mind-institute-biobank <?php echo ($expert_name == $ary['expert_name'] ? "current" : ""); ?>">
			<a href="index.php?expert_name=<?php echo $ary['expert_name']; ?>"><?php echo $ary['expert_name']; ?></a></li>
				<?php } ?>
		</ul>
		
        
	</li>
	
	<li class="healthy-brain-network">
		<a href="figures.php">Figures</a>
	</li>
	<!--
	<li class="healthy-brain-network">
		<a href="help.php">Help</a>
	</li>
	-->
	
	<li class="healthy-brain-network">
		<a href='https://cmi.hackpad.com/CMI-Annotator-PRbcvRqgaxj'>HackPad</a>
	</li>		

		
		</nav>

	</div>

	<div class="column primary">
		<section class="heading">
		<h2>
			
				
			
				Annotator
			
			
		</h2>
		</section>

<?php if (!$expert_name) { ?>
	Please choose an expert from the navigation menu on the left.
<?php } else { ?>
		<aside class="widgets">
			<div class="utils"></div>
	
							
<section class="module standalone parents-guide" style="text-align:center;margin: 0 auto;">
        
			<p class="anno_heading">
			<h3>Text Search</h3>
			<form id='freetextSelector' action='index.php' method='get'>
				<input type='hidden' name='expert_name' value='<?php echo $expert_name; ?>'>
				<?php if (isset($year)) { ?><input type='hidden' name='year' value='<?php echo $year; ?>'><?php } ?>
				<?php if (isset($tag)) { ?><input type='hidden' name='tag' value='<?php echo $tag; ?>'><?php } ?>
				<?php if (isset($status)) { ?><input type='hidden' name='status' value='<?php echo $status; ?>'><?php } ?>
				<?php if (isset($rsFmri)) { ?><input type='hidden' name='rsFmri' value='<?php echo $rsFmri; ?>'><?php } ?>
				<?php if (isset($clinicalTag)) { ?><input type='hidden' name='clinicalTag' value='<?php echo $clinicalTag; ?>'><?php } ?>				
				<input name='freetext' type='text' size="20" value="<?php if (isset($freetext)) echo $freetext; ?>"><br>
				<input name='freetextSubmit' type='submit' value='Go'>
			</form>		
			</p>	
			<p class="anno_heading">
			<h3>Year</h3>
			<form id='yearSelector' action='index.php' method='get'>
				<input type='hidden' name='expert_name' value='<?php echo $expert_name; ?>'>
				<?php if (isset($tag)) { ?><input type='hidden' name='tag' value='<?php echo $tag; ?>'><?php } ?>
				<?php if (isset($status)) { ?><input type='hidden' name='status' value='<?php echo $status; ?>'><?php } ?>
				<?php if (isset($freetext)) { ?><input type='hidden' name='freetext' value='<?php echo $freetext; ?>'><?php } ?>
				<?php if (isset($rsFmri)) { ?><input type='hidden' name='rsFmri' value='<?php echo $rsFmri; ?>'><?php } ?>
				<?php if (isset($clinicalTag)) { ?><input type='hidden' name='clinicalTag' value='<?php echo $clinicalTag; ?>'><?php } ?>													
				<select name='year' onchange='this.form.submit()'>
					<?php 
					echo "<option value=''";						 
					if (!isset($year) || "" == $year) { echo " selected"; } 
					echo ">(any)</option>"; 		
					for ($i = 1997; $i <= 2013; $i++) { 
						echo "<option value='" . $i . "'";						 
						if (isset($year) && $i == $year) { echo " selected"; } 
						echo ">" . $i . "</option>"; 			
					} 
					?>					
				</select>
			</form>
			</p>
			<p class="anno_heading">
			<h3>Tag</h3>
			<form id='tagSelector' action='index.php' method='get'>
				<input type='hidden' name='expert_name' value='<?php echo $expert_name; ?>'>
				<?php if (isset($year)) { ?><input type='hidden' name='year' value='<?php echo $year; ?>'><?php } ?>
				<?php if (isset($status)) { ?><input type='hidden' name='status' value='<?php echo $status; ?>'><?php } ?>
				<?php if (isset($freetext)) { ?><input type='hidden' name='freetext' value='<?php echo $freetext; ?>'><?php } ?>
				<?php if (isset($rsFmri)) { ?><input type='hidden' name='rsFmri' value='<?php echo $rsFmri; ?>'><?php } ?>					
				<select name='tag' onchange='this.form.submit()'>
					<?php 
					echo "<option value=''";						 
					if (!isset($tag) || "" == $tag) { echo " selected"; } 
					echo "></option>"; 	
					foreach ($tags as $tagId => $tagName) {
						echo "<option value='" . $tagId . "'";						 
						if (isset($tag) && $tagId == $tag) { echo " selected"; } 
						echo ">" . prettyTag($tagName) . "</option>"; 			
					} 
					?>					
				</select>
			</form>
			</p>
			<p class="anno_heading">
			<h3>Status</h3>
			<form id='statusSelector' action='index.php' method='get'>
				<input type='hidden' name='expert_name' value='<?php echo $expert_name; ?>'>
				<?php if (isset($year)) { ?><input type='hidden' name='year' value='<?php echo $year; ?>'><?php } ?>
				<?php if (isset($tag)) { ?><input type='hidden' name='tag' value='<?php echo $tag; ?>'><?php } ?>
				<?php if (isset($freetext)) { ?><input type='hidden' name='freetext' value='<?php echo $freetext; ?>'><?php } ?>
				<?php if (isset($rsFmri)) { ?><input type='hidden' name='rsFmri' value='<?php echo $rsFmri; ?>'><?php } ?>
				<?php if (isset($clinicalTag)) { ?><input type='hidden' name='clinicalTag' value='<?php echo $clinicalTag; ?>'><?php } ?>				
				<select name='status' onchange='this.form.submit()'>
					<?php 
					# Any
					echo "<option value='any'";						 
					if (!isset($status) || "any" == $status) { echo " selected"; } 
					echo ">(any)</option>";
					# Pending tagging					 
					echo "<option value='pending_tagging'";						 
					if (isset($status) && "pending_tagging" == $status) { echo " selected"; } 
					echo ">Pending Tagging</option>";
					# Tagging complete
					echo "<option value='tagging_complete'";						 
					if (isset($status) && "tagging_complete" == $status) { echo " selected"; } 
					echo ">Tagging Complete</option>";
					?>					
				</select>
			</form>		
			</p>
			<p class="anno_heading" style="text-align:center;">
			<h3>RS/fMRI</h3>
			<form id='rsFmriSelector' action='index.php' method='get'>
				<input type='hidden' name='expert_name' value='<?php echo $expert_name; ?>'>
				<?php if (isset($year)) { ?><input type='hidden' name='year' value='<?php echo $year; ?>'><?php } ?>
				<?php if (isset($tag)) { ?><input type='hidden' name='tag' value='<?php echo $tag; ?>'><?php } ?>
				<?php if (isset($freetext)) { ?><input type='hidden' name='freetext' value='<?php echo $freetext; ?>'><?php } ?>
				<?php if (isset($status)) { ?><input type='hidden' name='status' value='<?php echo $status; ?>'><?php } ?>
				<?php if (isset($clinicalTag)) { ?><input type='hidden' name='clinicalTag' value='<?php echo $clinicalTag; ?>'><?php } ?>				
				<label class="checkboxes" style='text-align:center;<?php if (isset($rsFmri)) { echo "font-weight: bold;"; } ?>'><input class="checkboxes" id="rsFmri" type='checkbox' name='rsFmri' onchange='this.form.submit()' <?php if (isset($rsFmri)) { echo "checked"; } ?> value='1'>&nbsp;<span class="checkboxes"><?php echo "Hide Results Tagged Not&nbsp;RS/Not&nbsp;fMRI"; ?></span></label>
			</form>		
			</p>				
			<?php if (count($breaks) > 0) { ?>
			<p class="anno_heading">
			<h3>Navigate</h3>			
			<form id='breakSelector' action='index.php' method='get'>
				<input type='hidden' name='expert_name' value='<?php echo $expert_name; ?>'>
				<?php if (isset($rsFmri)) { ?><input type='hidden' name='rsFmri' value='<?php echo $rsFmri; ?>'><?php } ?>
				<?php if (isset($year)) { ?><input type='hidden' name='year' value='<?php echo $year; ?>'><?php } ?>
				<?php if (isset($tag)) { ?><input type='hidden' name='tag' value='<?php echo $tag; ?>'><?php } ?>	
				<?php if (isset($status)) { ?><input type='hidden' name='status' value='<?php echo $status; ?>'><?php } ?>	
				<?php if (isset($freetext)) { ?><input type='hidden' name='freetext' value='<?php echo $freetext; ?>'><?php } ?>
				<?php if (isset($clinicalTag)) { ?><input type='hidden' name='clinicalTag' value='<?php echo $clinicalTag; ?>'><?php } ?>
				<select name='break' onchange='this.form.submit()'>
					<?php for ($i = 0; $i < count($breaks); $i++) {
						$text = ($breaks[$i]+1) . " - " . (($i+1) < count($breaks) ? ($breaks[$i+1]) : $overallCount);
						echo "<option value='" . $breaks[$i] . "'";
						if ($breaks[$i] == $break) {
							echo " selected";
						} 
						echo ">" . $text . "</option>"; 	
					} ?>					
				</select>
			</form>
			<?php } ?>
			</p>
			<p class="anno_heading">	
			<h3>Clinical Tag</h3>
			<form id='editClinicalTags' action='index.php' method='post'>
				<input type='hidden' name='expert_name' value='<?php echo $expert_name; ?>'>
				<?php if (isset($rsFmri)) { ?><input type='hidden' name='rsFmri' value='<?php echo $rsFmri; ?>'><?php } ?>
				<?php if (isset($year)) { ?><input type='hidden' name='year' value='<?php echo $year; ?>'><?php } ?>
				<?php if (isset($status)) { ?><input type='hidden' name='status' value='<?php echo $status; ?>'><?php } ?>
				<?php if (isset($freetext)) { ?><input type='hidden' name='freetext' value='<?php echo $freetext; ?>'><?php } ?>		
				<select name='clinicalTag' onchange='this.form.submit()'>
					<?php 
					echo "<option value=''";						 
					if (!isset($clinicalTag) || "" == $clinicalTag) { echo " selected"; } 
					echo "></option>"; 	
					foreach ($clinicalTags as $tagId => $tagName) {
						echo "<option value='" . $tagId . "'";						 
						if (isset($clinicalTag) && $tagId == $clinicalTag) { echo " selected"; } 
						echo ">" . prettyTag($tagName) . "</option>"; 			
					} 
					?>					
				</select><br>
				<input name='deleteClinicalTagSubmit' type='submit' value='Delete Selected'><br>
				<input name='addClinicalTag' type='text' size='20' value=''>&nbsp;<input name='addClinicalTagSubmit' type='submit' value='Add New'>
			</form>			


</section>

					
				
	
		</aside>

		<section class="module story">
 		<p align='center'>
 			<?php echo "Results " . min($break+1, $break+count($refs)) . " - " . ($break + count($refs)) . " of " . $overallCount . 
				(isset($tag) || isset($clinicalTag) ? ("<br> tagged " . $completeTagText) : "") . 				
				(isset($year) ? ("<br> in year " . $year) : "") . 
				(isset($status) && "any" != $status ? ("<br> with status " . prettyTag($status)) : "") .
				(isset($freetext) ? ("<br> matching text '" . $freetext . "'") : "") .
				(isset($rsFmri) ? ("<br> not tagged 'Not Resting State' or 'Not fMRI'") : ""); ?></b>
				</p>
		<?php 
		for ($i = 0; $i < count($refs); $i++) {
			$ref = $refs[$i];
			$nextIx = ($i == count($refs) - 1) ? ($i-1) : ($i); 
			$chosenTags = getChosenTagsString($ref['assignment_id']);
			?>		
		<div style='width:100%; border-top:1px solid #cccccc;'>			
				<a id='<?php echo "aid" . $i; ?>'></a>
				<table border=0 cellspacing=3 cellpadding=3>
				<?php foreach ($ref as $key => $value) { ?>
					<?php if (!in_array($key, $ignoreFields) && $value != "") { ?>
					<tr>
					<?php 
					$style = "";
					if ($key == "title" || $key == "author") {
						$style = "<b>";
					} 
					if ($key == "status") {
						$value = prettyTag($value);
					}
					if ($key == "pmid") {
						$value = "<a target='_blank' href='http://www.ncbi.nlm.nih.gov/pubmed/?term=" . $value . "%5Buid%5D'>" . $value . "</a>";
					}
					if ($key == "doi") {
						$value = "<form id='doi" . $i . "' target='_blank' method='post' action='http://dx.doi.org/'><input type='hidden' name='hdl' value='" . $value . "'> <a href='javascript:void(0);' onclick='document.getElementById(\"doi" . $i . "\").submit();'>" . $value . "</a></form>";
					}
					if ($key == "mesh_heading_list") {
						$key = "MeSH";
					}
					echo "<td width='15%' valign='top'>" . $style . prettyTag($key) . ":</td><td>" . $style . $value . "</td>"; ?>
					</tr>
					<?php } ?>
				<?php } ?>	
				</table>
				<table><tr><td>	
						<?php echo getCheckForm($nextIx); ?>
						<input type='hidden' name='ref_id' value='<?php echo $ref['ref_id']; ?>'>
						<input type='hidden' name='expert_name' value='<?php echo $expert_name; ?>'>
						<input type='hidden' name='assignment_id' value='<?php echo $ref['assignment_id']; ?>'>
						<?php if (isset($year)) { ?><input type='hidden' name='year' value='<?php echo $year; ?>'><?php } ?>				
						<?php if (isset($tag)) { ?><input type='hidden' name='tag' value='<?php echo $tag; ?>'><?php } ?>	
						<?php if (isset($status)) { ?><input type='hidden' name='status' value='<?php echo $status; ?>'><?php } ?>	
						<?php if (isset($freetext)) { ?><input type='hidden' name='freetext' value='<?php echo $freetext; ?>'><?php } ?>
						<?php if (isset($clinicalTag)) { ?><input type='hidden' name='clinicalTag' value='<?php echo $clinicalTag; ?>'><?php } ?>
						<?php if (isset($rsFmri)) { ?><input type='hidden' name='rsFmri' value='<?php echo $rsFmri; ?>'><?php } ?>									

					<div class='checkboxes' style='float:left; width:100%; text-align:left;'>
						<table style='border-spacing:0px; border-collapse: separate;'>
						<tr><td style='padding:0px;vertical-align:middle'>
						Age Range:</td><td style='padding:3px;vertical-align:middle'><input id='ageLow' size='5' type='text' name='ageRangeLow' value='<?php echo $ref['age_range_low']; ?>'> to <input id='ageHigh' size='5' type='text' name='ageRangeHigh' value='<?php echo $ref['age_range_high']; ?>'>
						</td></tr>
						<tr><td style='padding:0px;vertical-align:middle'>
						Sample Size:</td><td style='padding:3px;vertical-align:middle'><input id='sampleSize' size='5' type='text' name='sampleSize' value='<?php echo $ref['sample_size']; ?>'>
						</td></tr>
						</table> 
					</div>
					<div class='checkboxes' style='float:left; width:50%; text-align:left; margin-top:5px;'>
						<p class="anno_heading">
						Issues
						<?php echoTagHtml(55, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(56, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(61, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(97, $tags, $chosenTags, 2); ?>
						</p>						
						<p class="anno_heading">
						General Tags					
						<?php echoTagHtml(25, $tags, $chosenTags, 2); ?> 	
						<?php echoTagHtml(26, $tags, $chosenTags, 4); ?>
						<?php echoTagHtml(27, $tags, $chosenTags, 4); ?>
						<?php echoTagHtml(7, $tags, $chosenTags, 4); ?>
						<?php echoTagHtml(9, $tags, $chosenTags, 4); ?>
						<?php echoTagHtml(17, $tags, $chosenTags, 4); ?>
						<?php echoTagHtml(71, $tags, $chosenTags, 4); ?>
						<?php echoTagHtml(46, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(47, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(48, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(49, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(60, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(90, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(91, $tags, $chosenTags, 2); ?>
						</p>
						<p class="anno_heading">
						Analysis Methods
						<?php echoTagHtml(50, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(51, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(52, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(53, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(54, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(84, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(85, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(86, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(95, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(96, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(77, $tags, $chosenTags, 2); ?>
						</p>
						<p class="anno_heading">
						Other Imaging Modalities
						<?php echoTagHtml(45, $tags, $chosenTags, 2); ?>
						<?php echoTagHtml(57, $tags, $chosenTags, 4); ?>
						<?php echoTagHtml(58, $tags, $chosenTags, 4); ?>
						<?php echoTagHtml(59, $tags, $chosenTags, 4); ?>						
						<?php echoTagHtml(62, $tags, $chosenTags, 4); ?>
						<?php echoTagHtml(63, $tags, $chosenTags, 4); ?>
						<?php echoTagHtml(81, $tags, $chosenTags, 4); ?>												
						<?php echoTagHtml(92, $tags, $chosenTags, 4); ?>
						<?php echoTagHtml(93, $tags, $chosenTags, 4); ?>
						<?php echoTagHtml(94, $tags, $chosenTags, 4); ?>
						<?php echoTagHtml(76, $tags, $chosenTags, 4); ?>
						</p>						
					</div>
					<div class='checkboxes' style='float:left; width:50%; text-align:left; margin-top:5px;'>
						<p>
						<?php echoTagHtml(20, $tags, $chosenTags); ?>		
						<?php foreach ($clinicalTags as $tagId => $tagName) {							
							echoTagHtml($tagId, $clinicalTags, $chosenTags, 2);  
						} ?>
						</p>
					</div>					
					<div style='float:left; width:90%; text-align:center;'>
					<?php if ($ref['status'] == 'tagging_complete') { ?>
						<input type='submit' name='mark_refresh_button' value='Mark PENDING TAGGING'>
					<?php } else if ($ref['status'] == 'pending_tagging') { ?>	
						<input type='submit' name='mark_refresh_button' value='Mark TAGGING COMPLETE'>
					<?php } ?>
					</div>						
			</form>
			</td></tr></table>
		</div>
		<p></p>
		<?php } ?>	

 	
		</section>
<?php } ?>
		
		
	</div>



			</div>
		</div> <!--! end of #wrapper -->

		<footer>
			<a href="http://www.childmind.org/en/healthy-brain-network/cmi-librarian/#modal_guide" id="guide-link" class="modal" style="display: none;"></a>
			<div class="wrapper">
				<nav>
                    
					<ul>
						<li><a href="http://www.childmind.org/">Home</a></li>
						
						<li><a href="http://www.childmind.org/press/">Press Room</a></li>
						
						<li><a href="http://www.childmind.org/contact/">Contact Us</a></li>
						
						<li><a href="http://www.childmind.org/privacy/">Privacy</a></li>
						
						<li><a href="http://www.childmind.org/site-map/">Site Map</a></li>
						
						<li><a href="http://www.childmind.org/terms-of-use/">Terms of Use</a></li>
						
					</ul>
                    
					<small>Copyright 2014 All rights reserved
</small>
					<ul>
						<li><strong>General Inquiries Call 212.308.3118</strong></li>
						<li>445 Park Avenue, New York, NY 10022</li>
						<li>
							<ul class="social-icons">
								<li><a href="http://feeds.feedburner.com/ChildMindDotOrg" target="_blank"><img src="./cmi_librarian_files/rss.png" alt="RSS" width="20" height="20"></a></li>
								<li><a href="http://www.facebook.com/ChildMindInstitute" target="_blank"><img src="./cmi_librarian_files/facebook.png" alt="Facebook" width="20" height="20"></a></li>
								<li><a href="http://twitter.com/ChildMindDotOrg" target="_blank"><img src="./cmi_librarian_files/twitter.png" alt="Twitter" width="20" height="20"></a></li>
								<li><a href="https://plus.google.com/u/0/b/104916968679971130579/" target="_blank"><img src="./cmi_librarian_files/gplus.png" alt="Google+" width="20" height="20"></a></li>
								<li><a href="http://www.flickr.com/people/54841878@N07/" target="_blank"><img src="./cmi_librarian_files/flickr.png" alt="Flickr" width="20" height="20"></a></li>
								<li><a href="http://www.youtube.com/childmindinstitute" target="_blank"><img src="./cmi_librarian_files/youtube.png" alt="YouTube" width="20" height="20"></a></li>
							</ul>
						</li>
					</ul>
				</nav>
			</div>
		</footer>

		<div id="load-spinner"></div>

		<div id="background-shim" class="science-innovation"></div>


</body></html>