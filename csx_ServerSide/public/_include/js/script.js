$(document).ready(function() {

    $("#fit_range").prop('disabled', true);

    // form validation


    $('#fit_box').change(function() {
        if (this.checked) {
            $("#fit_range").prop('disabled', false);
        } else {
            $("#fit_range").prop('disabled', true);
        }
    });

});

