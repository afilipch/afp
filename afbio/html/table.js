function sortTable(n, data_type) {
  var table, rows, switching, i, x, y, shouldSwitch, dir, switchcount = 0;
  table = document.getElementById('myTable');
  switching = true;
  //Set the sorting direction to ascending:
  dir = 'asc'; 
  /*Make a loop that will continue until
  no switching has been done:*/
  while (switching) {
    //start by saying: no switching is done:
    switching = false;
    rows = table.rows;
    /*Loop through all table rows (except the
    first, which contains table headers):*/
    for (i = 1; i < (rows.length - 1); i++) {
      //start by saying there should be no switching:
      shouldSwitch = false;
      /*Get the two elements you want to compare,
      one from current row and one from the next:*/
      x = rows[i].getElementsByTagName('TD')[n];
      y = rows[i + 1].getElementsByTagName('TD')[n];
      /*check if the two rows should switch place,
      based on the direction, asc or desc:*/
      if (dir == 'asc') {
        if (data_type == 0 && x.innerHTML.toLowerCase() > y.innerHTML.toLowerCase()) {
          //if so, mark as a switch and break the loop:
          shouldSwitch = true;
          break;
        }
        else if (data_type == 1 && Number(x.innerHTML) > Number(y.innerHTML)) {
            shouldSwitch = true;
            break;
        }
      } else if (dir == 'desc') {
        if (data_type == 0 && x.innerHTML.toLowerCase() < y.innerHTML.toLowerCase()) {
          //if so, mark as a switch and break the loop:
          shouldSwitch = true;
          break;
        }
        else if (data_type == 1 && Number(x.innerHTML) < Number(y.innerHTML)) {
            shouldSwitch = true;
            break;
        }
      }
    }
    if (shouldSwitch) {
      /*If a switch has been marked, make the switch
      and mark that a switch has been done:*/
      rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
      switching = true;
      //Each time a switch is done, increase this count by 1:
      switchcount ++;      
    } else {
      /*If no switching has been done AND the direction is 'asc',
      set the direction to 'desc' and run the while loop again.*/
      if (switchcount == 0 && dir == 'asc') {
        dir = 'desc';
        switching = true;
      }
    }
  }
}


function my_search(n, input_name, table_name) {
  // Declare variables
  var input, filter, table, tr, td, i, txtValue;
  input = document.getElementById(input_name);
  filter = input.value.toUpperCase();
  table = document.getElementById(table_name);
  tr = table.getElementsByTagName("tr");

  // Loop through all table rows, and hide those who don't match the search query
  for (i = 0; i < tr.length; i++) {
    td = tr[i].getElementsByTagName("td")[n];
    if (td) {
      txtValue = td.textContent || td.innerText;
      if (txtValue.toUpperCase().indexOf(filter) > -1) {
        tr[i].style.display = "";
      } else {
        tr[i].style.display = "none";
      }
    }
  }
}


function my_filter_greater(n, input_name, table_name) {
  // Declare variables
  var input, filter, table, tr, td, i, txtValue;
  input = document.getElementById(input_name);
  filter = Number(input.value)
  table = document.getElementById(table_name);
  tr = table.getElementsByTagName("tr");

  // Loop through all table rows, and hide those who don't match the search query
  for (i = 0; i < tr.length; i++) {
    td = tr[i].getElementsByTagName("td")[n];
    if (td) {
      txtValue = td.textContent || td.innerText;
      if (Number(txtValue) > filter) {
        tr[i].style.display = "";
      } else {
        tr[i].style.display = "none";
      }
    }
  }
}


function my_filter_lesser_abs(n, input_name, table_name) {
  // Declare variables
  var input, filter, table, tr, td, i, txtValue;
  input = document.getElementById(input_name);
  filter = Number(input.value)
  table = document.getElementById(table_name);
  tr = table.getElementsByTagName("tr");

  // Loop through all table rows, and hide those who don't match the search query
  for (i = 0; i < tr.length; i++) {
    td = tr[i].getElementsByTagName("td")[n];
    if (td) {
      txtValue = td.textContent || td.innerText;
      if (Math.abs(Number(txtValue)) < filter) {
        tr[i].style.display = "";
      } else {
        tr[i].style.display = "none";
      }
    }
  }
}



