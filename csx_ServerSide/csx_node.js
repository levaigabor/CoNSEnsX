
// Node.js goodies <3

var express    = require('express');
var fs         = require("fs");
var bodyParser = require('body-parser');
var multer     = require('multer');

// creating express application
var app = express();

app.use(express.static('public'));
app.use(bodyParser.urlencoded({ extended: false }));
app.use(multer({ dest: 'uploads/'}));

// serving CSX start page on GET
app.get('/index.htm', function (req, res) {
   res.sendFile( __dirname + "/" + "index.html" );
})


app.post('/CSX_start', function (req, res) {
  // ID generation for CSX
  function makeid() {
    var text = "";
    var possible = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";

    for( var i=0; i < 6; i++ )
        text += possible.charAt(Math.floor(Math.random() * possible.length));

    return text;
  }

  var calculation_ID = makeid();

  // CSX console command, will be make up according to the posted form
  var csx_call = "consensx.py -i " + calculation_ID + " ";

  var upLoadFile = function(file){
    fs.readFile(file.path, function(error, data) {
      fs.writeFile(file, data, function(error) {
        if (error) {
          console.log(error);
        }
      });
    });
  }

  for (var i = 0; i < req.files.length; i++) {
    upLoadFile(req.files[i]);
  }

  // PDB ID if given
  var pdb_web = req.body.pdb_web;

  if (pdb_web != undefined) {
    // if PDB ID is given, append it to the CSX command
    csx_call += "-p " + pdb_web + " ";
  } else {
    // if PDB file is given, save it and append it to the CSX command
    csx_call += "-f uploads/" + req.files.pdb_upload.name + " ";
  }

  // append STAR-NMR file to the CSX command
  csx_call += "-b " + req.files.bmrb_upload.name + " ";

  var NOE_upload = req.files.NOE_upload

  // if PDB file is given, save it and append it to the CSX command
  if (NOE_upload != undefined) {
    csx_call += "-r " + req.files.NOE_upload.name + " ";
  }

  // append KARPLUS parameter set selection to the CSX command
  csx_call += "-d " + req.body.KARPLUS + " ";

  // if given, add SVD option
  var RDCSVD = req.body.RDCSVD;

  if (RDCSVD != undefined) {
    csx_call += "-R ";
  }

  // RDC LC model
  csx_call += "-l " + req.body.RDCLC + " ";

  // if given, add r3 averaging option
  var r3average = req.body.r3average;

  if (r3average != undefined) {
    csx_call += "-r3 ";
  }

  var superimpose = req.body.superimpose;
  var fit_range = req.body.fit_range;

  if (superimpose != undefined) {
    csx_call += "-s ";
    if (fit_range != undefined) {
      csx_call += "--fit_range " + fit_range + " ";
    }
  }

  console.log(csx_call);


  res.end( JSON.stringify(  ) );
})





var server = app.listen(8081, function () {

  var host = server.address().address
  var port = server.address().port

  console.log("Example app listening at http://%s:%s", host, port)

})
