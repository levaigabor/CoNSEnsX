$(document).ready(function() {
	$( "#selection").click(function () {
		$(".involve-result").toggleClass ( "involve", 1000, "easeOutSine")
	});
	
	$( ".inv-1" ).click(function() {
  		$(".group3").toggleClass( "involve", 1000, "easeOutSine" );
	});

	$( ".inv-2" ).click(function() {
  		$(".group4").toggleClass( "involve", 1000, "easeOutSine" );
	});

	$( ".inv-3" ).click(function() {
  		$(".group5").toggleClass( "involve", 1000, "easeOutSine" );
	});

	$( ".inv-4" ).click(function() {
  		$(".group6").toggleClass( "involve", 1000, "easeOutSine" );
	});

});
