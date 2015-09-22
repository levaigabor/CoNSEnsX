$(document).ready(function() {
	$( "#selection").click(function () {
		$(".involve-result").toggleClass ( "involve", 1000, "easeOutSine")
	});
	
	$( ".show" ).click(function() {
  		$(".group3").toggleClass( "involve", 1000, "easeOutSine" );
	});

});
