(function () {
  "use strict";

  document.addEventListener("DOMContentLoaded", function () {
    var selects = document.querySelectorAll(".blend-doc-select");
    var i;

    for (i = 0; i < selects.length; i++) {
      selects[i].addEventListener("change", function () {
        var base;
        var candidate;
        var currentBase = null;
        var j;
        var relativePath;
        var target;

        if (!this.value) {
          return;
        }

        for (j = 0; j < this.options.length; j++) {
          candidate = new URL(this.options[j].value, window.location.href);
          if (candidate.origin === window.location.origin &&
              window.location.pathname.indexOf(candidate.pathname) === 0 &&
              (!currentBase || candidate.pathname.length > currentBase.pathname.length)) {
            currentBase = candidate;
          }
        }

        target = new URL(this.value, window.location.href);
        if (currentBase) {
          base = currentBase.pathname;
          relativePath = window.location.pathname.substring(base.length);
          target = new URL(relativePath, target);
          target.search = window.location.search;
          target.hash = window.location.hash;
        }

        window.location.href = target.href;
      });
    }
  });
}());
