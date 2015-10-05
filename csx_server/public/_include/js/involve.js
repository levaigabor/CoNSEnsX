$(document).ready(function() {
    $("#selection").click(function () {
        $(".involve-result").toggleClass ( "involve", 1000, "easeOutSine")
    });
});

// button naming
// add radio buttons
// check before upload!!

function myFunction() {
    var ranges = document.querySelectorAll(".inputrange");
    var measures = document.querySelectorAll(".options");
    var measure;
    var range_selected = false;
    var command = "";
    var i;

    // add selected measure to command
    for (i = 0; i < measures.length; i++) {
        if (measures[i].checked) {
            console.log("selected measure: " + measures[i].value);
            command += "MEASURE " + measures[i].value + "\n";
            measure = measures[i].value;
        }
    }

    // check if measure is selected
    if (measure == undefined) {
        alert("Please select compliance measure!");
        return;
    }

    // add involved parameters to command
    for (i = 0; i < ranges.length; i++) {
        if (ranges[i].value != 0) {
            command += ranges[i].id + " " + ranges[i].value + "\n";
            range_selected = true;
        }
    }

    // check if at least one parameter is involved in selection
    if (!range_selected) {
        alert("Please invole at least one parameter!");
        return;
    }

    // UPLAD HERE!
    console.log(command);
}
