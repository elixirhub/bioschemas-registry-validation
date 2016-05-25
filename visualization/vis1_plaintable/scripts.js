/**
 * Created by robertopreste on 25/05/2016.
 */

function tableMaker(csvfile) {
    // Create a table reporting all the scraped websites and their types

    var table = d3.select("#table1")
        .append("table");
    
    var thead = table.append("thead");
    var tbody = table.append("tbody");
    var columns = ["Website", "Event", "Organization", "Person", "Training"];
    
    thead.append("tr")
        .selectAll("th")
        .data(columns)
        .enter()
        .append("th")
        .text(function (columns) {
            return columns;
        });

    d3.csv(csvfile, function (error, data) {
        data = data.map(function (d) {
            d.value = [d["Website"], d["Event"], d["Organization"], d["Person"], d["Training"]];
            return d;
        });

        var rows = tbody.selectAll("tr")
            .data(data)
            .enter()
            .append("tr");

        rows.selectAll("td")
            .data(function (row) {
                return columns.map(function (column) {
                    return {column: column, value: row[column]};
                });
            })
            .enter()
            .append("td")
            .html(function (d) {
                return d.value;
            });
    });

}

tableMaker("test1.csv");