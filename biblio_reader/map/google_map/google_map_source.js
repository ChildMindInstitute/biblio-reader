function initMap() {
        var uluru = {lat: -25.363, lng: 131.044};
        var map = new google.maps.Map(document.getElementById('map'), {
          zoom: 2,
          center: {lat: 0, lng: 0}
        });
        var marker = new google.maps.Marker({
          position: uluru,
          map: map
          })
        var marker = new google.maps.Marker({
          position: {lat: 42.3631542, lng:-71.0688334},
          map: map
          })
        $getJSON()
}