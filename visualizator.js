/**
 * Created by robertopreste on 21/06/2016.
 */

// Main page

function main_page() {
    var tot_props = document.getElementById("tot_props");
    var tot_web = document.getElementById("tot_websites");
    var tot_type = document.getElementById("tot_types");
    var last_up = document.getElementById("last_update");

    d3.csv("props_reg.csv", function(error, data) {
        tot_props.innerHTML = data.length;
    });

    d3.csv("websites_reg.csv", function(error, data) {
        tot_web.innerHTML = data.length;
    });

    d3.csv("types_reg.csv", function(error, data) {
        tot_type.innerHTML = data.length;
    });

    var today = new Date();
    var dd = today.getDate();
    var mm = today.getMonth() + 1;
    var yyyy = today.getFullYear();
    if (dd < 10) { dd = '0'+dd; }
    if (mm < 10) { mm = '0'+mm; }

    today = mm + '/' + dd + '/' + yyyy;
    last_up.innerHTML = today;
}



// Types Registry Page
/**
var tableBody = document.getElementById("types_table");
        d3.json("../registry.json", function(data) {
            for (var i in Object.keys(data)) {
                //var row = tableBody.insertRow(0);
                var el = Object.keys(data)[i];
                console.log(el);    // Prima cella (Website)
                //console.log(typeof data[el]);
                for (var j in Object.keys(data[el])) {
                    var il = Object.keys(data[el][j]);
                    for (var k = 0; k < il.length; k++) {
                        console.log(il[k]); // Seconda cella (Type)
                        console.log(data[il[k]]);
                        //for (var l = 0; l < il[k].length; l++) {
                        //    console.log(il[k][l]);
                        //}
                    }
                }

            }

        });
*/
/**
var tableBody = document.getElementById("types_table");
d3.json("../registry.json", function(data) {
    for (var i in Object.keys(data)) {
        var el = Object.keys(data)[i];
        console.log(el);
        console.log(typeof data[el]);
        for (var j in Object.keys(data[el])) {
            var il = Object.keys(data[el][j]);
            for (var k = 0; k < il.length; k++) {
                console.log(il[k]);
            }
        }

    }

});
*/

/**
var table = d3.select("#types_table").append("table");
var thead = table.append("thead");
var tbody = table.append("tbody");
var columns = ["Type", "Websites", "Entries"];

var listtest = [];
var listtest2 = [];

thead.append("tr")
    .selectAll("th")
    .data(columns)
    .enter()
    .append("th")
    .text(function(columns) { return columns; });

d3.csv("../types_reg.csv", function(error, data) {
    for (var i of data) {
        console.log(i.Types);
        if (i == "Person") {
            console.log("hello person");
        }
    }
});

console.log(listtest2);
console.log(listtest2.length);
console.log(listtest);
console.log(listtest.length);


var types_table = document.getElementById("types_table");
var tbody = types_table.getElementsByTagName("tbody");
var row = tbody.insertRow(0);
var cell1 = row.insertCell(0);
var cell2 = row.insertCell(1);

d3.csv("types_reg.csv", function(error, data) {
    var newText = document.createTextNode(data[1]);
    cell1.appendChild(newText);
});
*/
