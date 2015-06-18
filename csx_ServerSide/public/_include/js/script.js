$(document).ready(function() {

    $("#fit_range").prop('disabled', true);

    // form validation
    $('#submit_button').click(function() {
        var hasPDB  = ($('#pdb_file').val() != "" ) || ($('#pdb_id').val() != "" );
        var hasMBMR = $('#bmbr_file').val() != ""

        if ((hasPDB == false) || (hasMBMR == false)) {
            alert("You have to select a PDB and a BMBR file at least to start calculation");
            return false;
        }

        var hasBothPDB = ($('#pdb_file').val() != "" ) && ($('#pdb_id').val() != "" );

        if (hasBothPDB) {
            alert("Please deselect one of the PDF inputs");
        }

        return false;
    });

    $('#fit_box').change(function() {
        if (this.checked) {
            $("#fit_range").prop('disabled', false);
        } else {
            $("#fit_range").prop('disabled', true);
        }
    });

});
