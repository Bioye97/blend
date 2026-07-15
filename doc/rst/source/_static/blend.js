(function () {
  "use strict";

  document.addEventListener("DOMContentLoaded", function () {
    var selects = document.querySelectorAll(".blend-doc-select");
    var i;

    for (i = 0; i < selects.length; i++) {
      selects[i].addEventListener("change", function () {
        if (this.value) {
          window.location.href = this.value;
        }
      });
    }
  });
}());
