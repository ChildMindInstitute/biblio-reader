/*
File: global.js

About: Version
	1.0

Project: Child Mind Institute

Description:
	A common file that includes all globally shared functionality for CMI

Requires:
	- jQuery <http://jquery.com/>

*/

/*
Class: CMI
	Scoped to the CMI Global Namespace
*/
var CMI = CMI || {};

// When the DOM is ready.
(function () {

	// Storing a variable to reference
	var $self = CMI;

	/*
	Namespace: CMI.vars
		Shared global variables
	*/
	$self.vars = {

		/*
		variable: queue
			Contains the functions ready to be fired on DOM ready
		*/
		queue : []
	};

	/*
	Namespace: CMI.legacy
		A legacy namespace to be used if CMI has any legacy dependencies
	*/
	$self.legacy = {};

	/*
	Namespace: CMI.utils
		Shared global utilities
	*/
	$self.utils = {
		/*
		sub: ie6Check
			Adds an IE6 flag to jQuery
		*/
		ieCheck : function() {
			// Let's set a flag for IE 6
			$.extend($.browser, {
				ie6 : function () {
					return !!($.browser.msie && $.browser.version == 6);
				}()
			});
		}(),

		/*
		sub: queue
			A global initializer. Takes a function argument and queues it until <init> is fired

		Parameters:
			object - The function to initialize when the DOM is ready
		*/
		queue : function (object) {
			$self.vars.queue.push(object);
		},

		/*
		sub: init
			When fired, loops through $self.vars.queue and fires each queued function
		*/
		init : function() {
			var queue = $self.vars.queue;

			$.each(queue, function(i, object) {
				for (var key in object) {
					if (object.hasOwnProperty(key) && (typeof object[key] === "function")) {
						object[key]();
					}
				}
			});
		},

		// cookie helper functions
		createCookie : function (name, value, days) {
			if (days) {
				var date = new Date();
				date.setTime(date.getTime()+(days*24*60*60*1000));
				var expires = "; expires="+date.toGMTString();
			}
			else var expires = "";
			document.cookie = name+"="+value+expires+"; path=/";
		},

		eraseCookie : function (name) {
			$self.utils.createCookie(name,"",-1);
		},

		readCookie : function (name) {
			var nameEQ = name + "=";
			var ca = document.cookie.split(';');
			for(var i=0;i < ca.length;i++) {
				var c = ca[i];
				while (c.charAt(0)==' ') c = c.substring(1,c.length);
				if (c.indexOf(nameEQ) == 0) return c.substring(nameEQ.length,c.length);
			}
			return null;
		},

		/**
		* Shamelessly ripped from Modernizr <http://modernizr.com>
		*
		* test_props is a generic CSS / DOM property test; if a browser supports
		*   a certain property, it won't return undefined for it.
		*   A supported CSS property returns empty string when its not yet set.
		*/
		testProps : function( props, callback ) {
			var mod = 'modernizr',
					doc = document,
					m = doc.createElement( mod ),
					m_style = m.style;

			for ( var i in props ) {
				if ( m_style[ props[i] ] !== undefined && ( !callback || callback( props[i], m ) ) ) {
					return true;
				}
			}

			return false;
		},

			testPropsAll : function( prop, callback ) {
			var uc_prop = prop.charAt(0).toUpperCase() + prop.substr(1),

			// following spec is to expose vendor-specific style properties as:
			//	 elem.style.WebkitBorderRadius
			// and the following would be incorrect:
			//	 elem.style.webkitBorderRadius
			// Webkit and Mozilla are nice enough to ghost their properties in the lowercase
			//	 version but Opera does not.

			// see more here: http://github.com/Modernizr/Modernizr/issues/issue/21
			props = [
				prop,
				'Webkit' + uc_prop,
				'Moz' + uc_prop,
				'O' + uc_prop,
				'ms' + uc_prop,
				'Khtml' + uc_prop
			];

			return !!this.testProps( props, callback );
		},

		testTransitions : function() {
			return this.testPropsAll("transitionProperty");
		},

		test3d : function() {
			//  set_css_all( 'perspective:500' );

			var ret = !!this.testProps([ 'perspectiveProperty', 'WebkitPerspective', 'MozPerspective', 'OPerspective', 'msPerspective' ]),
					prefixes = ' -o- -moz- -ms- -webkit- -khtml- '.split(' ');

			// webkit has 3d transforms disabled for chrome, though
			//   it works fine in safari on leopard and snow leopard
			// as a result, it 'recognizes' the syntax and throws a false positive
			// thus we must do a more thorough check:
			if (ret){
				var st = document.createElement('style'),
						doc = document,
						docElement = doc.documentElement,
						div = doc.createElement('div');

				// webkit allows this media query to succeed only if the feature is enabled.
				// "@media (transform-3d),(-o-transform-3d),(-moz-transform-3d),(-ms-transform-3d),(-webkit-transform-3d),(modernizr){#modernizr{height:3px}}"
				st.textContent = '@media ('+prefixes.join('transform-3d),(')+'modernizr){#modernizr{height:3px}}';
				doc.getElementsByTagName('head')[0].appendChild(st);
				div.id = 'modernizr';
				docElement.appendChild(div);

				ret = div.offsetHeight === 3;

				st.parentNode.removeChild(st);
				div.parentNode.removeChild(div);
			}
			return ret;
		},

		PageController : function(control) {

			var Controller = function(control) {
				return this.init(control);
			};

			Controller.prototype = {
				_timer : 1000,

				generateController : function(control, list, items) {
					var _self = this;

					// Shell list item
					var ul = $('<ul></ul>');

					// Generate an icon for each list item
					items.each(function() {
						var item = $('<li><a href="#"></a></li>');
						item.attr("class", $(this).attr("class"));

						item.bind("click", function(e) {
							e.preventDefault();

							var icons = ul.children().not(".page-control-previous").not(".page-control-next"),
									idx = icons.index(this);

							$(".page-control-active", control).removeClass("page-control-active").find("h3").animate({
								opacity : 0
							}, {
								duration : _self._timer
							});

							$(this).addClass("page-control-active");
							items.eq(idx).addClass("page-control-active").find("h3").animate({
								opacity : 1
							}, {
								duration : _self._timer
							});

							_self.positionControllerBasedOn(list, items, items.eq(idx));
						});

						ul.append(item);
					});

					ul.find("li").bind("click", function(e, triggered) {
						if (!triggered) {
							_self.autoStop();
						}
					});

					// Prevent Link Defaults
					control.find("li").bind("click", function(e, triggered) {

						var el = $(this);

						if (!el.is(".page-control-playing")) {

							if (!triggered) {
								_self.autoStop();
							}

							var yt = $(this).attr("data-unique-id");
							if (yt) {
								_self.playInlineVideo(control, $(this), yt);
							}
						}
					});

					control.find("a[href*='youtube.com']").bind("click", function(e) {
						e.preventDefault();
					});

					var triggerNext = function(type) {
						var icons = ul.children().not(".page-control-previous").not(".page-control-next"),
								active = icons.filter(".page-control-active"),
								next = active[type]();

						if (!next.attr("class")) {
							next.trigger("click");
						}
					};

					$('<li class="page-control-previous"></li>').bind("click", function() {
						triggerNext("prev");
					}).prependTo(ul);

					$('<li class="page-control-next"></li>').bind("click", function() {
						triggerNext("next");
					}).appendTo(ul);

					ul.find("a").bind("click", function(e) {
						e.preventDefault();
					});

					// Append the list to the controller wrapper
					// Append the controller wrapper to control
					$('<div class="page-control-controller"></div>').append(ul).appendTo(control);

					var totalWidth = 0;

					ul.find("li").each(function(i, item) {
						totalWidth += $(this).outerWidth(true);
					});

					ul.parent().width(totalWidth);
				},

				positionControllerBasedOn: function(list, items, active) {
					var index = items.index(active);

					list.animate({
						left : -(active.outerWidth(true) * index),
						leaveTransforms : true
					}, {
						duration : this._timer
					});
				},

				autoStart : function(control) {
					var items = $(".page-control-controller li", control).not(".page-control-previous").not(".page-control-previous"),
							next;

					this._interval = window.setInterval(function() {
						next = items.filter(".page-control-active").next("li:not(.page-control-previous):not(.page-control-next)");

						if (!next.length) {
							next = items.eq(0);
						}

						next.trigger("click", true);
					}, 5000);
				},

				autoStop : function() {
					if (this._interval) {
						window.clearInterval(this._interval);
					}
				},

				playInlineVideo : function(control, active, id) {
					var steps, i, j,
							url = "http://www.youtube.com/v/" + id + "?fs=1&amp;hl=en_US",
							controller = control.find(".page-control-controller"),
							nodes = active.find("img, .page-control-lower-third"),
							_interval = 350;

					function activateVideo() {
						var link = url + '?fs=1&amp;hl=en_US&amp;rel=0&amp;autoplay=1';
						if (!active.find(".page-control-video").get(0)) {
							active.append('<div class="page-control-video"><object width="480" height="300"><param name="movie" value="' + link + '" /><param name="allowFullScreen" value="true" /><param name="wmode" value="transparent" /><param name="allowscriptaccess" value="always" /><embed src="' + link + '" type="application/x-shockwave-flash" allowscriptaccess="always" allowfullscreen="true" wmode="transparent" width="480" height="300" /></object><div class="btn-close-video"></div></div>');
						}
					}

					active.addClass("page-control-playing");

					if (!$self.utils.test3d()) {
						active.addClass("page-control-play");

						nodes.animate({
							left : 664
						}, {
							duration : 500,
							complete : activateVideo
						});
					} else {
						active.addClass("page-control-3d");

						steps = [
							function() {
								active.addClass("page-control-scale");
							},

							function() {
								active.removeClass("page-control-scale");
								active.addClass("page-control-rotate");
							},

							function() {
								active.removeClass("page-control-rotate");
								active.addClass("page-control-play");
							}
						];

						for (i = 0, j = steps.length; i < j; i++) {
							window.setTimeout(steps[i], _interval * i);
						}

						activateVideo();
					}

					controller.animate({
						top : 100
					}, {
						duration : 500
					});

					active.find(".btn-close-video").live("click", function() {
						if (!$self.utils.test3d()) {
							nodes.animate({
								left : 0
							}, {
								duration : 500,
								complete : function() {
									active.find(".page-control-video").remove();
									active.removeClass("page-control-playing");
									active.removeClass("page-control-play");
								}
							});
						} else {
							steps = [
								function() {
									active.removeClass("page-control-play");
									active.addClass("page-control-rotate");
								},

								function() {
									active.removeClass("page-control-rotate");
									active.addClass("page-control-scale");
								},

								function() {
									active.removeClass("page-control-scale");
								},

								function() {
									active.find(".page-control-video").remove();
									active.removeClass("page-control-playing");
								}
							];

							for (i = 0, j = steps.length; i < j; i++) {
								window.setTimeout(steps[i], _interval * i);
							}
						}

						controller.animate({
							top : 0
						}, {
							duration : 500
						});
					});
				},

				init : function(control) {
					// Too hot for iOS
					if ("ontouchstart" in window) {
						return;
					}

					var _self = this;

					var content = control.find(".page-control-content"),
							list = content.find("> ul"),
							items = list.children(),
							active = items.filter(".page-control-active").eq(0) || items.eq(0);

					if (!active.get(0)) {
						active = items.eq(0).addClass("page-control-active");
					}

					// Sanity check. Escape out if less than two feature items
					if (items.length <= 1) {
						return;
					}

					// Set list width wide enough to have all items in a row.
					list.width(items.eq(0).outerWidth(true) * items.length);

					// Hide visibility on titles
					items.not(active).find("h3").css("opacity", 0);

					// Generate paginator
					_self.generateController(control, list, items);

					// Initial positioning
					_self.positionControllerBasedOn(list, items, active);

					// Auto Start
					if (items.length > 2) {
						_self.autoStart(control);
					}
				}
			};

			return Controller;
		}()
	};

	/*
	Namespace: CMI
		Under the CMI Local Namespace
	*/

	/*
	Function: global
		Takes care of a few global functionalities.
	*/
	$self.global = function () {
		if ($.browser.msie) {
			try {
				// Enable background image cache
					document.execCommand("BackgroundImageCache", false, true);
			} catch (ex) {}
		}

		// Set default easing
		$.easing.def = "easeInOutExpo";

		// JS is enabled
		$("html").removeClass("no-js");

		$("input[type='text'][placeholder]").each(function() {
			if (!("placeholder" in this)) {
				var input = $(this),
						placeholder = input.attr("placeholder");

				if (!input.val()) {
					input.val(placeholder);
				}

				input.bind({
					focus : function() {
						if (input.val() === placeholder) {
							input.val("");
						}
					},

					blur : function() {
						if (!input.val()) {
							input.val(placeholder);
						}
					}
				});
			}
		});

		// Unfix the background if less than 578px
		$(window).bind("resize", function() {
			if ($(window).height() <= 578) {
				$("body").addClass("unfix");
			} else {
				$("body").removeClass("unfix");
			}
		}).trigger("resize");
	};

	/*
	Function: header
		Encapsulates functionality found in the #header space
	*/
	$self.header = function() {
		$("header select").bind("change", function() {
			window.location.href = $(this).val();
		});

		var signUp = $("#email-sign-up"),
				panel = $("#donation-panel"),
				overlay = $("#email-overlay"),
				oWidth, oLeft, _timer = 500,
				reset = overlay.find("input[type='reset']"),
				oForm = overlay.find("form");

		signUp.bind("click", function(e) {
			e.preventDefault();

			panel.addClass("active");
			oWidth = oWidth || $(this).outerWidth();
				oLeft = oLeft || $(this).position().left;

			var inputTest = document.createElement("input");
			if (!("placeholder" in inputTest)) {
				$("input[type='text']", overlay).each(function() {
					var input = $(this),
							placeholder = input.attr("placeholder");

					if (!input.val()) {
						input.val(placeholder);
					}

					input.bind({
						focus : function() {
							if (input.val() === placeholder) {
								input.val("");
							}
						},

						blur : function() {
							if (!input.val()) {
								input.val(placeholder);
							}
						}
					});
				});
			}

			overlay.css({
				width : oWidth,
				left : oLeft
			}).fadeOut(0).fadeIn(_timer * 0.25, function() {
				$(this).animate({
					left : 0,
					width : "100%"
				}).find("form").fadeIn(_timer, function() {
					$(this).find("input").eq(0).focus();
				});
			}).find("form").fadeOut(0);
		});

		reset.bind("click", function(e) {
			e.preventDefault();

			oForm.fadeOut(_timer);

			overlay.animate({
				width : oWidth,
				left : oLeft
			}, {
				duration : _timer,
				complete : function() {
					overlay.fadeOut(_timer * 0.25, function() {
						panel.removeClass("active");
					});
				}
			});
		});

		// Chrome on windows can't yet handle inset shadows with rounded corners
		// No way to feature detect that, so we'll unfortunately have to
		// Disable inset shadows for Chrome/Win
		// See: media/resources/sass/src/modues/_header.scss
		if ((/win/i).test(window.navigator.platform) && (/chrome/i).test(window.navigator.userAgent)) {
			$("header #donation-panel input[type='text']").addClass("chrome-win");
		}

		/*
		Prompt users running old IE versions to upgrade… to yet another shitty IE release. ¯\_( ツ)_/¯
		Rather than the preferred method of using an SSI in the template, we have to include this here
		so we can cookie a dismissal in order to suppress the warning on subsequent page views.
		*/
		if ($.browser.msie && $.browser.version < 9) {
			var cookies = $self.utils.readCookie("CMI-IE_Update");
			if (cookies != "declined") {
				$("#wrapper").prepend("<div id='luddites'><a href='http://windows.microsoft.com/en-US/internet-explorer/products/ie/home' target='_blank'>You are using an outdated version of Internet Explorer. For a faster, safer browsing experience upgrade today.</a><a href='#close'>Close</a></div>");

				$("#luddites a[href='#close']").live("click", function () {
					$self.utils.createCookie("CMI-IE_Update", "declined", 1);
					$("#luddites").remove();
				});
			}
		}

	};

	/*
	Function: actionable
		Handles .actionable-button functionality
	*/
	$self.actionable = function() {
		var buttons = $(".actionable-button"),
				button, select, defaultText;

		if ($.browser.msie && $.browser.version < 8) {
			return;
		}

		buttons.each(function(i) {
			button = $(this);

			if (button.is("select")) {
				var wrap = $('<div></div>');
				wrap.attr("class", button.attr("class")).addClass("select");

				if (button.attr("disabled")) {
					wrap.addClass("disabled");
				}

				wrap.css("max-width", button.width());

				button.removeClass("actionable-button");
				button.wrap(wrap);

				if ((button.attr("selectedIndex") === 0) && !button.find(":selected").get(0).getAttribute("selected") && button.find("optgroup").get(0)) {
					defaultText = button.find("optgroup").attr("label");
				} else {
					defaultText = button.find(":selected").text();
				}

				$('<span>' + defaultText.toUpperCase() + '</span>').insertBefore(button);

				button.bind("change", function() {
					var el = $(this);
					el.prev().text(el.find(":selected").text().toUpperCase());
				});
			}
		});
	};

	/*
	Function: modules
		Handles module functionality
	*/
	$self.modules = function() {
		$("section").each(function () {
			if ($(this).children().html() == null) {
				$(this).css("display", "none");
			}
		});

		$(".module select").each(function() {
			var el = $(this),
					node = ($.browser.msie && $.browser.version < 8) ? el : el.parent(),
					submit = node.siblings(":submit");

			if (submit.get(0)) {
				el.bind("change", function() {
					$(this).closest("form").submit();
				});

				submit.remove();
			}
		});

		$(".alt-modules").each(function () {
			var moduleHeight = 0, headingHeight = 0,
					items = $(this).find("> li");

			items.each(function() {
				headingHeight = Math.max(headingHeight, $("> h3, > strong", this).height());
			}).find("> h3, > strong").css("min-height", headingHeight + 1);

			items.each(function() {
				moduleHeight = Math.max(moduleHeight, $(this).height());
			}).css("min-height", moduleHeight + 1);
		});

		$(".module.features, .module.standalone.related").each(function() {
			var moduleHeight = 0, headingHeight = 0,
					items = $(this).find("> ul > li");

			items.each(function() {
				headingHeight = Math.max(headingHeight, $("> h3, > strong", this).height());
			}).find("> h3, > strong").css("min-height", headingHeight + 1);

			items.each(function() {
				moduleHeight = Math.max(moduleHeight, $(this).height());
			}).css("min-height", moduleHeight + 1);
		});

		$(".threeCols, .module.body, .related ul, #widgets").each(function () {
			var moduleHeight = 0, headingHeight = 0, equalHeight = 0,
				items = $(this).children(".module");

			items.each(function() {
				headingHeight = Math.max(headingHeight, $("> h3, > strong", this).height());
			}).find("> h3, > strong").css("min-height", headingHeight + 1);

			items.find(".equalize").each(function() {
				equalHeight = Math.max(equalHeight, $(this).height());
			}).css("min-height", equalHeight + 1);

			items.each(function() {
				moduleHeight = Math.max(moduleHeight, $(this).height());
			}).css("min-height", moduleHeight + 1);
		});
	};

	/*
	Function: pageControl
		Global carousel module functionality
	*/
	$self.pageControl = function() {
		$(".page-control").each(function() {
			var control = new $self.utils.PageController($(this));
		});
	};

	/*
	Function: accordionControl
		Global accordion module functionality
	*/
	$self.accordionControl = function() {
		var controls = $(".accordion-control"),
				_timer = 500;

		controls.each(function() {
			var control = $(this),
					list = $("> dl", control),
					toggles = $("> dt", list),
					items = $("> dd", list);

			control.addClass("enabled");

			toggles.bind("click", function() {
				var el = $(this),
						active = el.hasClass("active"),
						item = el.next("dd");

				toggles.removeClass("active");
				items.removeClass("active");

				if (!active) {
					el.addClass("active");
					item.addClass("active");

					item.slideDown(_timer);
					items.not(item).slideUp(_timer);
				} else {
					item.slideUp(_timer);
				}
			});

			items.slideUp(0);

			window.setTimeout(function() {
				toggles.eq(0).trigger("click");
			}, _timer);
		});

		$(".expand-all").bind("click", function(e) {
			e.preventDefault();
			$(".accordion-control dd").slideDown(0);
		});
	};

	/*
	Function: tabControl
		Global tab control module functionality
	*/
	$self.tabControl = function() {
		$(".tab-control").each(function () {
			var el = $(this),
					controller = el.find(".tab-control-controller"),
					content = el.find(".tab-control-content");

			controller.children().bind("click", function (e) {
				var relative = content.find("." + $(this).attr("class"));
				relative.show().siblings().hide();

				$(this).addClass("active").siblings().removeClass("active");
			}).find("a").bind("click", function (e) {
				e.preventDefault();
			});

			controller.children().eq(0).trigger("click");

			content.css("min-height", controller.height());
		});
	};

	$self.isLoading = (function () {
		var $body = $("body");

		return function (isLoading) {
			if (isLoading) {
				$body.addClass("loading");
			} else {
				$body.removeClass("loading");
			}
		};
	}());

	/*
	Function: modalControl
		Global modal content control
	*/

	$self.modalControl = (function() {

		var positionModal = function (container) {
			container = container || $("#lightbox");

			container.css("top", function () {
				var winHeight = $(window).height(),
					modalHeight = container.outerHeight();

				return Math.max(Math.floor((winHeight - modalHeight) / 2), 20);
			}).show();
		};

		var openModal = function (contentEl, modalClass) {
			var content, darkroom, lightbox;
			modalClass = modalClass || "non";

			// Chrome only supports inset shadows with rounded corners on Mac OS X
			// No way to feature detect that, so unfortunately we'll have to disable inset shadows for Chrome Win/Lin
			// See: media/resources/sass/src/modues/_modals.scss
			if (!(/mac/i).test(window.navigator.platform) && (/chrome/i).test(window.navigator.userAgent)) {
				$("#modals input").addClass("chrome-fix");
			}

			content = $(contentEl).html();
			darkroom = $("<div id='darkroom'></div>");
			lightbox = $("<div id='lightbox' class='lightbox " + modalClass +"'></div>");

			lightbox.fadeOut().append(content);
			darkroom.append(lightbox);
			$("body").append(darkroom);

			darkroom.fadeIn(400);
			positionModal(lightbox);

			$(document).keyup(function (e) {
				var keyCode = e.keyCode || window.event.keyCode; // window.event is for IE
				if (keyCode == 27) { // esc key closes modal
					closeModal();
				}
			});
		};

		var closeModal = function () {
			$("#darkroom").fadeOut(400, function () {
				$(this).remove();
			});

			if (!$(this).hasClass("sufk-submit") || !$(this).hasClass("bts-submit")) {
				return false;
			}
		};

		$("a.modal").bind("click", function (e) {
			var selector = $(this).attr("href");
			openModal(selector, selector.replace("#modal_", ""));
			e.preventDefault();
		});

		$("#lightbox .close, .sufk-submit").live("click", closeModal);

		return {
			position: positionModal,
			open: openModal,
			close: closeModal
		};
	}());

	$self.showGuide = function () {
		var $link = $("a#guide-link"),
		$form = $("#form-box > form");

		if ($link.length !== 0 && !$.cookie("download-guide")) {
			$link.trigger("click");
			// set for 2 weeks by default
			$.cookie("download-guide", "yes", {expires: 14});

			$form.live("submit", function (e) {
				$.cookie("download-guide", "yes", { expires: 365 });
				$("#darkroom").fadeOut(400, function () {
					$(this).remove();
				});
			});
		}
	};

	$self.checkerSurveys = function () {
		var path = window.location.pathname;
		if (path.search("symptom-checker") !== -1 &&
			path.search("results") === -1) {
			return;
		}

		var leftKey = "left-symptom-checker";
		var completedKey = "completed-symptom-checker";
		var leftEarly = $.cookie(leftKey);
		var completedChecker = $.cookie(completedKey);
		var isResults = $("body").hasClass("page-results");

		var bindModalForm = function (selector) {
			var $modal = $(".lightbox");
			var $form = $(selector);
			var $content = $modal.find(".content");

			$form.submit(function (e) {
				var inputData = $form.serializeArray();
				var data = {};

				$.each(inputData, function (i, input) {
					if (input.value) {
						data[input.name] = input.value;
					}
				});

				$form.find(":input").attr("disabled", true);
				$self.isLoading(true);

				$.ajax({
					url: $form.attr("action"),
					type: $form.attr("method"),
					cache: false,
					data: data,
					contentType: "application/json; charset=utf-8",
					dataType: "json",
					success: function (response) {
						$modal.find(".close").text("Close");
						$content.html("<p>" + response.msg + "</p>");
						$self.modalControl.position();
						$self.isLoading(false);
					}
				});

				return false;
			});
		};

		if (leftEarly || completedChecker || isResults) {
			$(".survey-alert").show().bind("click", function () {
				var modal = leftEarly ? "#modal_survey-incomplete" : "#modal_survey-evaluation";
				var form = leftEarly ? ".survey-form-incomplete" : ".survey-form-complete";

				$(this).hide();
				$self.modalControl.open(modal, "survey");
				bindModalForm(form);

				$.removeCookie(leftKey, {path: "/"});
				$.removeCookie(completedKey, {path: "/"});
			});
		}
	};

	/*
	Function: printControl
		Global print hijax control
	*/
	$self.printControl = function() {
		// Enable printer links
		$("li.print > a, a.print").bind("click", function () {
			window.print();
			return false;
		});
	};

	/*
	Function: filterControl
		Global filter control
	*/
	$self.filterControl = function () {
		var control = $("#filter"),
				selects = control.find("select");

		selects.bind("change", function () {
			control.trigger("submit");
		});
	};

	/*
	Function: reportLinks
		Global comment reporting links control
	*/
	$self.reportLinks = function() {
		$("a.report").bind("click", function () {
			var link = $(this);

			$.ajax({
				type: "POST",
				url: "/report/"+ $(link).attr("rel"),
				success: function () {
					$(link).removeClass("report").addClass("reported").html("This comment has been reported");
				}
			});
			return false;
		});
	};

	$self.externalLinks = function () {
		var host = window.location.hostname.split(".").slice(1).join(".");
		$("a[href^='http']:not([href*='" + host + "'])").each(function () {
			$(this).attr("target", "_blank");
		});
	};

	$self.charCounter = function () {
		var max = $("#id_text").attr("maxlength");
		$("p.char_counter").html(max + ' characters remaining');

		$("#id_text").keyup(function () {
			var rem,str = $(this).val();

			if (max - str.length <= 5) {
				rem = "<em>"+ (max - str.length) +"</em>";
			} else {
				rem = max - str.length;
			}

			$("p.char_counter").html(rem + ' characters remaining');
		});
	};

	/*
	Callback: queue
		Sends local functions to a global queuer for initialization See: <CMI.utils.queue>
	*/
	$self.utils.queue($self);

}).call(CMI);