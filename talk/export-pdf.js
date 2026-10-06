const { chromium } = require('playwright');
(async () => {
  const b = await chromium.launch({executablePath: '/opt/pw-browsers/chromium'});
  const p = await b.newPage({viewport: {width: 1280, height: 720}});
  await p.goto('file://' + __dirname + '/slides.html?print-pdf', {waitUntil: 'networkidle'});
  await p.waitForTimeout(3000);
  await p.pdf({path: '' + __dirname + '/slides.pdf', width: '1280px', height: '720px', printBackground: true});
  await b.close();
})();
