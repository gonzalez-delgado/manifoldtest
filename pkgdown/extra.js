// Collapsible function families on the reference index page.
//
// pkgdown renders each `subtitle` of _pkgdown.yml as a section holding the
// heading and its description, followed by a section holding the list of
// functions. Each pair is folded into a <details> block, collapsed by default,
// whose header row is the arrow followed by the family name.
//
// The styles live here rather than in extra.css so that the markup and its
// layout always ship together (a cached stylesheet cannot break the layout).
var REF_FAMILY_CSS = [
  ".ref-family-toggle { margin-bottom: 1rem; }",
  "details.ref-family {",
  "  border: 1px solid var(--bs-border-color, #dee2e6);",
  "  border-radius: 8px;",
  "  margin-bottom: 0.75rem;",
  "}",
  "details.ref-family > summary {",
  "  display: block;",
  "  list-style: none;",
  "  cursor: pointer;",
  "  padding: 0.75rem 1rem;",
  "  border-radius: 8px;",
  "}",
  "details.ref-family > summary::-webkit-details-marker { display: none; }",
  "details.ref-family > summary::marker { content: none; }",
  "details.ref-family > summary:hover { background-color: var(--bs-tertiary-bg, #f8f9fa); }",
  ".ref-family-title { display: flex; align-items: center; gap: 0.6rem; }",
  ".ref-family-title h3 { margin: 0; font-size: 1.25rem; }",
  ".ref-family-arrow {",
  "  flex: none;",
  "  width: 0.5rem;",
  "  height: 0.5rem;",
  "  border-right: 2px solid currentColor;",
  "  border-bottom: 2px solid currentColor;",
  "  transform: rotate(-45deg);",
  "  transition: transform 0.15s ease;",
  "}",
  "details.ref-family[open] .ref-family-arrow { transform: rotate(45deg); }",
  "details.ref-family > summary .section-desc p { margin: 0.25rem 0 0 1.1rem; color: #6c757d; }",
  "details.ref-family > dl {",
  "  margin: 0;",
  "  padding: 0.5rem 1rem 0.75rem 2.1rem;",
  "  border-top: 1px solid var(--bs-border-color, #dee2e6);",
  "}"
].join("\n");

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

    // Header row: arrow and family name on one line. The heading keeps its
    // id, so links and the "On this page" menu still point to it.
    var summary = document.createElement("summary");
    var title = document.createElement("div");
    title.className = "ref-family-title";
    var arrow = document.createElement("span");
    arrow.className = "ref-family-arrow";
    arrow.setAttribute("aria-hidden", "true");
    h3.querySelectorAll("a.anchor").forEach(function (a) { a.remove(); });
    title.appendChild(arrow);
    title.appendChild(h3);
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

  var style = document.createElement("style");
  style.textContent = REF_FAMILY_CSS;
  document.head.appendChild(style);

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
  var pageHeader = main.querySelector(".page-header");
  pageHeader.parentNode.insertBefore(toggle, pageHeader.nextSibling);

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
