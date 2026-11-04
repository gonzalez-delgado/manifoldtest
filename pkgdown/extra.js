// Collapsible function families on the reference index page.
//
// pkgdown renders each `subtitle` of _pkgdown.yml as a section holding the
// heading and its description, followed by a section holding the list of
// functions. Each pair is folded into a <details> block, collapsed by default.
document.addEventListener("DOMContentLoaded", function () {
  if (!document.querySelector(".template-reference-index")) return;
  var main = document.querySelector("main");
  if (!main) return;

  var families = [];
  main.querySelectorAll("div.section.level2 > h3").forEach(function (h3) {
    var header = h3.parentElement;
    var body = header.nextElementSibling;
    var list = body && body.querySelector("dl");
    if (!list || list.children.length === 0) return;

    var details = document.createElement("details");
    details.className = "ref-family";

    // The count sits next to the heading, not inside it, so that the
    // "On this page" menu built from the headings shows the title only
    var summary = document.createElement("summary");
    var title = document.createElement("div");
    title.className = "ref-family-title";
    var count = list.querySelectorAll("dt").length;
    var badge = document.createElement("span");
    badge.className = "ref-family-count";
    badge.textContent = count + (count === 1 ? " function" : " functions");
    h3.querySelectorAll("a.anchor").forEach(function (a) { a.remove(); });
    title.appendChild(h3);
    title.appendChild(badge);
    summary.appendChild(title);
    var desc = header.querySelector(".section-desc");
    if (desc) summary.appendChild(desc);

    details.appendChild(summary);
    details.appendChild(list);
    header.querySelectorAll("dl").forEach(function (dl) { dl.remove(); });
    header.appendChild(details);
    body.remove();
    families.push(details);
  });
  if (families.length === 0) return;

  // Expand / collapse all
  var toggle = document.createElement("button");
  toggle.type = "button";
  toggle.className = "btn btn-sm btn-outline-secondary ref-family-toggle";
  var sync = function () {
    var allOpen = families.every(function (d) { return d.open; });
    toggle.textContent = allOpen ? "Collapse all" : "Expand all";
  };
  toggle.addEventListener("click", function () {
    var open = !families.every(function (d) { return d.open; });
    families.forEach(function (d) { d.open = open; });
    sync();
  });
  families.forEach(function (d) { d.addEventListener("toggle", sync); });
  sync();
  var header = main.querySelector(".page-header");
  header.parentNode.insertBefore(toggle, header.nextSibling);

  // Open the family targeted by a link (e.g. from the "On this page" menu)
  var openTarget = function () {
    var id = decodeURIComponent(window.location.hash.slice(1));
    var target = id && document.getElementById(id);
    var family = target && target.closest("details.ref-family");
    if (family) {
      family.open = true;
      family.scrollIntoView();
    }
  };
  window.addEventListener("hashchange", openTarget);
  openTarget();
});
