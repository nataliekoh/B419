$(document).ready(function() {
   'use strict';

    var format = function(d) {
        return d3.format('.03%')(d);
    };

    var prev_years = [2005, 2010];

    prev_years.forEach(function(year) {
        console.log(year);
        var map = d3.geomap.choropleth()
            .geofile('d3-geomap/topojson/world/countries.json')
            .column('Prevalence (' + year + ')')
            .format(format)
            .legend(true)
            .unitId('country code');

        d3.csv('data/tb_prev_' + year + '.csv', function(error, data) {
            d3.select('#prev_map_' + year)
                .datum(data)
                .call(map.draw, map);
        });
    });
});