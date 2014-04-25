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
$ignoreFields = array('ref_id', 'tags', 'assignment_id', 'master_tags', 'expert_name', 'assigned', 'author', 'affiliation', 'age_range_low', 'age_range_high', 'update_time', 'ta_id', 'ta_assignments_id', 'ta_tags_id', 'volume', 'pages', 'number', 'isbn', 'journal_iso');			
$maxResults = 50;

mysql_connect('localhost', 'root', 'renew');
mysql_select_db('annotate');

$expert_name = "";

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
	
	<li class="healthy-brain-network current" style="background-color:white">
		<a href="index.php" style="background-color:white">Annotator</a>
        
        
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
	
	<li class="healthy-brain-network current">
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
			
				
			
				Figures
			
			
		</h2>
		</section>

		<section class="module story" style="text-align:center;width:710px">

			<h1>
			<p><a name="piles">Annotation Status</a></p>
			<p><a href="../figures/piles.pdf"><img src="../figures/piles.png" width=710></a></p>				
			<p><a name="overallgrowth">Overall Growth</a></p>
			<p><a href="../figures/overall_growth.pdf"><img src="../figures/overall_growth.png"></a></p>		
			<p><a name="journaldist">Journal Distribution</a></p>
			<p><a href="../figures/journal_dist.pdf"><img src="../figures/journal_dist.png" width=710></a></p>		
			<p><a name="clinicaldist">Tag Distribution</a></p>
			<p><a href="../figures/clinical_bytag_hist.pdf"><img src="../figures/clinical_bytag_hist.png" width=710></a></p>		
			<p><a name="methodsgrowth">Methods Growth</a></p>
			<p><a href="../figures/methods_growth_hist.pdf"><img src="../figures/methods_growth_hist.png" width=710></a></p>
			<p><a name="tfidf">DF vs. Conditional TF</a></p>
			<p><a href="../figures/tfidf_top.pdf"><img src="../figures/tfidf_top.png" width=710></a></p>
			</h1>
			<p><h1>References</h1></p>
			<p>Doherty et al. "Bibliometric Analysis of Resting State."<br> 3rd Biennial Conference on Resting State Brain Connectivity, 2012. <a href="../figures/doherty_2012_conrs.pdf">[pdf]</a></p>
			<p>&nbsp;</p>			
		</section>
		
		
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